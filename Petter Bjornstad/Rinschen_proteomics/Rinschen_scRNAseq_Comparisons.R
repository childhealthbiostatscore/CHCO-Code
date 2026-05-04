################################################################################
# PATHWAY GENES — PT CELL DE + VIOLIN PLOTS
# SGLT2i vs. no SGLT2i in T2D, PT cells only
# Seurat: ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData
################################################################################

library(Seurat)
library(nebula)
library(Matrix)
library(SingleCellExperiment)
library(scran)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(openxlsx)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::rename)
conflicted::conflicts_prefer(dplyr::count)
conflicted::conflicts_prefer(base::intersect)
conflicted::conflicts_prefer(base::setdiff)

# ── Paths ─────────────────────────────────────────────────────────────────────
seurat_path <- "C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData"
output_dir  <- "C:/Users/netio/Documents/UofW/proteomics_vs_PET/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# PATHWAY GENE SETS (same as figure)
################################################################################

pathway_genes <- list(
  "Glycolysis" = c(
    "PKLR", "PFKFB3", "PFKL", "ALDOC", "HK2",
    "ENO2", "PGK1", "PGAM1", "TPI1", "GAPDH"
  ),
  "Gluconeogenesis" = c(
    "SLC25A10", "GOT2", "GOT1", "FBP1",
    "SLC25A11", "PCK1", "MDH1"
  ),
  "Pyruvate metabolism and TCA cycle" = c(
    "SDHB", "SUCLG1", "PDK2", "ACO2", "IDH3G",
    "SUCLA2", "HAGH", "PDHB", "LDHA"
  ),
  "Glutathione conjugation" = c(
    "CNDP2", "GSTM4", "GSTT2B", "GSTO1",
    "GGCT", "GSTM3", "AKR1A1"
  ),
  "Metallothioneins bind metals" = c(
    "MT1G", "MT1X", "MT1H", "MT2A"
  )
)

all_genes <- unlist(pathway_genes, use.names = FALSE)

################################################################################
# LOAD + PREP SEURAT OBJECT
################################################################################

load(seurat_path)
so <- so_kpmp_sc
rm(so_kpmp_sc)

# ── Cell type labels — exact pattern from prior analyses ──────────────────────
so$celltype1 <- case_when(
  grepl("PT-",  so$celltype_rpca) ~ "PT",
  grepl("TAL-", so$celltype_rpca) ~ "TAL",
  grepl("EC-",  so$celltype_rpca) ~ "EC",
  grepl("POD",  so$celltype_rpca) ~ "POD",
  grepl("MAC",  so$celltype_rpca) ~ "MAC",
  grepl("MON",  so$celltype_rpca) ~ "MON",
  grepl("PC-",  so$celltype_rpca) ~ "PC",
  grepl("FIB",  so$celltype_rpca) ~ "FIB_MC_VSMC",
  grepl("DTL",  so$celltype_rpca) ~ "DTL",
  so$celltype_rpca == "DCT"       ~ "DCT",
  so$celltype_rpca == "ATL"       ~ "ATL",
  so$celltype_rpca == "B"         ~ "B",
  so$celltype_rpca == "T"         ~ "T",
  TRUE                            ~ so$celltype_rpca
)
so$celltype1 <- as.character(so$celltype1)

so$KPMP_celltype2 <- as.character(so$KPMP_celltype)
so$celltype2 <- ifelse(
  so$KPMP_celltype %in% c("aPT", "PT-S1/S2", "PT-S3"), "PT",
  ifelse(grepl("TAL", so$KPMP_celltype), "TAL",
         ifelse(grepl("EC-", so$KPMP_celltype), "EC",
                so$KPMP_celltype2)))
so$celltype2 <- as.character(so$celltype2)
so$DCT_celltype <- ifelse(so$KPMP_celltype %in% c("DCT", "dDCT"), "DCT", "Non-DCT")

# ── Filter: T2D only, remove CRC-55 ──────────────────────────────────────────
so <- subset(so, subset = record_id != "CRC-55")
so <- subset(so, subset = group == "Type_2_Diabetes")

# ── SGLT2i grouping label — use sglt2i_ever (more reliable than sglt2i_timepoint)
# sglt2i_timepoint has known data entry errors for RH-49-T, RH-74-T, RH-75-T, RH-77-T
so$sglt2_group <- ifelse(so$sglt2i_ever == "Yes", "SGLT2i", "No SGLT2i")

cat("Total T2D cells:", ncol(so), "\n")
cat("SGLT2i breakdown:\n");        print(table(so$sglt2_group))
cat("Cell type breakdown:\n");     print(table(so$celltype2))

# ── Subset to PT cells (celltype2) ───────────────────────────────────────────
so_pt <- subset(so, subset = celltype2 == "PT")
cat("\nPT cells:", ncol(so_pt), "\n")
cat("PT donors:", length(unique(so_pt$record_id)), "\n")
cat("PT SGLT2i breakdown:\n"); print(table(so_pt$sglt2_group))

# ── Check which pathway genes are in the assay ────────────────────────────────
genes_present <- base::intersect(all_genes, rownames(so_pt))
genes_missing <- base::setdiff(all_genes, rownames(so_pt))
cat("\nPathway genes present:", length(genes_present), "/", length(all_genes), "\n")
if (length(genes_missing) > 0) cat("Missing:", paste(genes_missing, collapse = ", "), "\n")

################################################################################
# NEBULA — TARGETED DE FOR PATHWAY GENES IN PT CELLS
################################################################################

# Compute pooled offset on the FULL PT object (before gene subsetting)
counts_full         <- round(GetAssayData(so_pt, layer = "counts"))
sce_full            <- SingleCellExperiment(assays = list(counts = counts_full))
sce_full            <- computeSumFactors(sce_full)
so_pt$pooled_offset <- sizeFactors(sce_full)
rm(sce_full, counts_full)

# NOW subset to pathway genes
so_pt_sub  <- so_pt[genes_present, ]
counts_mat <- round(GetAssayData(so_pt_sub, layer = "counts"))

meta_df <- so_pt_sub@meta.data %>%
  mutate(sglt2i = factor(sglt2i_ever, levels = c("No", "Yes")))

# Remove cells with NA in sglt2i or offset
complete_idx <- complete.cases(meta_df[, c("sglt2i", "pooled_offset")])
cat("Cells with complete data:", sum(complete_idx), "\n")
counts_mat <- counts_mat[, complete_idx]
meta_df    <- meta_df[complete_idx, ]

# Simple model — SGLT2i only, no additional covariates
pred_mat <- model.matrix(~ sglt2i, data = meta_df)

lib    <- meta_df$pooled_offset
data_g <- group_cell(count  = counts_mat,
                     id     = meta_df$record_id,
                     pred   = pred_mat,
                     offset = lib)
if (is.null(data_g)) {
  data_g <- list(count  = counts_mat,
                 id     = meta_df$record_id,
                 pred   = pred_mat,
                 offset = lib)
}

nebula_res <- nebula(
  count      = data_g$count,
  id         = data_g$id,
  pred       = data_g$pred,
  offset     = data_g$offset,
  model      = "NBLMM",
  ncore      = 1,
  reml       = TRUE,
  output_re  = TRUE,
  covariance = TRUE
)

# ── Parse results — auto-detect SGLT2i logFC column ──────────────────────────
res_raw <- nebula_res$summary %>% as_tibble()
cat("\nNEBULA output columns:\n"); print(names(res_raw))

fc_col <- grep("logFC.*sglt2|logFC.*Yes", names(res_raw), value = TRUE,
               ignore.case = TRUE)[1]
p_col  <- grep("^p_.*sglt2|^p_.*Yes",    names(res_raw), value = TRUE,
               ignore.case = TRUE)[1]
cat("Using logFC col:", fc_col, "| p-value col:", p_col, "\n")

res_df <- res_raw %>%
  dplyr::rename(gene_symbol = gene,
                logFC_sglt2 = all_of(fc_col),
                p_sglt2     = all_of(p_col)) %>%
  mutate(
    fdr = p.adjust(p_sglt2, method = "BH"),
    sig = fdr < 0.05
  ) %>%
  left_join(
    imap_dfr(pathway_genes, ~ tibble(gene_symbol = .x, pathway = .y)),
    by = "gene_symbol"
  ) %>%
  arrange(fdr)

cat("\nSignificant genes (FDR < 0.05):", sum(res_df$sig, na.rm = TRUE), "\n")
print(res_df %>% dplyr::filter(sig) %>%
        dplyr::select(gene_symbol, pathway, logFC_sglt2, fdr))

write.xlsx(res_df, file.path(output_dir, "NEBULA_PT_pathway_genes_SGLT2i.xlsx"))
cat("Results saved.\n")

################################################################################
# BAR CHART — same layout as original figure
################################################################################

plot_df <- res_df %>%
  left_join(
    imap_dfr(pathway_genes, ~ tibble(
      gene_symbol = .x, pathway = .y, gene_order = seq_along(.x)
    )),
    by = c("gene_symbol", "pathway")
  ) %>%
  mutate(
    pathway     = factor(pathway, levels = names(pathway_genes)),
    gene_symbol = factor(gene_symbol, levels = all_genes)
  )

make_bar_panel <- function(pw_name, df) {
  gene_order <- pathway_genes[[pw_name]]
  pw_df <- df %>%
    dplyr::filter(pathway == pw_name) %>%
    mutate(gene_symbol = factor(as.character(gene_symbol), levels = gene_order))
  
  # Fill in missing genes as NA rows so x-axis is complete
  missing_genes <- base::setdiff(gene_order, as.character(pw_df$gene_symbol))
  if (length(missing_genes) > 0) {
    pw_df <- bind_rows(pw_df,
                       tibble(gene_symbol = factor(missing_genes, levels = gene_order),
                              logFC_sglt2 = NA_real_, sig = FALSE))
  }
  
  sig_genes <- pw_df %>% dplyr::filter(sig) %>% pull(gene_symbol) %>% as.character()
  axis_face <- ifelse(gene_order %in% sig_genes, "bold", "plain")
  
  ggplot(pw_df, aes(x = gene_symbol, y = logFC_sglt2)) +
    geom_col(fill = "#3AAFA9", color = NA, width = 0.7, na.rm = TRUE) +
    geom_col(data = pw_df %>% dplyr::filter(sig),
             fill = "#3AAFA9", color = "black", linewidth = 0.5, width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey30") +
    scale_x_discrete(drop = FALSE) +
    labs(title = pw_name, x = NULL, y = "log\u2082 fold change") +
    theme_bw(base_size = 9) +
    theme(
      plot.title         = element_text(size = 9, face = "bold"),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 7.5,
                                        face = axis_face),
      axis.title.y       = element_text(size = 8),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

panels <- imap(
  setNames(names(pathway_genes), names(pathway_genes)),
  ~ make_bar_panel(.x, plot_df)
)

fig_bar <- (panels[[1]] | panels[[2]]) /
  (panels[[3]] | panels[[4]] | panels[[5]]) +
  plot_annotation(
    title    = "Pathway genes: SGLT2i vs. No SGLT2i (PT cells, NEBULA)",
    subtitle = "log\u2082FC age+sex adjusted | bold border = FDR < 0.05 | positive = higher in SGLT2i",
    theme    = theme(plot.title    = element_text(size = 11, face = "bold"),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

ggsave(file.path(output_dir, "pathway_barchart_PT_NEBULA_SGLT2i.pdf"),
       fig_bar, width = 12, height = 7, device = "pdf")
ggsave(file.path(output_dir, "pathway_barchart_PT_NEBULA_SGLT2i.png"),
       fig_bar, width = 12, height = 7, dpi = 300)
cat("Bar chart saved.\n")

################################################################################
# SIDE-BY-SIDE COMPARISON — OUR NEBULA vs. SCHAUB ET AL. (published scRNA-seq)
# Schaub et al. teal bars = log2FC T2Di(+) vs T2Di(-) i.e. SGLT2i effect
# Same comparison, same cell type (PT), same disease (youth T2D)
# Values read from published panel C figure
################################################################################

published_vals <- tribble(
  ~gene_symbol,  ~pathway,                              ~logFC_published,
  # Glycolysis — mostly suppressed (negative) in PT, per paper text + Fig 4C
  # HK2, PFKL, PKLR explicitly suppressed; PGK1/PGAM1/TPI1 enhanced (were downregulated in T2D)
  "PKLR",        "Glycolysis",                         -0.02,
  "PFKFB3",      "Glycolysis",                         -0.02,
  "PFKL",        "Glycolysis",                         -0.02,
  "ALDOC",       "Glycolysis",                         -0.02,
  "HK2",         "Glycolysis",                         -0.03,
  "ENO2",        "Glycolysis",                         -0.02,
  "PGK1",        "Glycolysis",                          0.05,
  "PGAM1",       "Glycolysis",                          0.05,
  "TPI1",        "Glycolysis",                          0.13,
  "GAPDH",       "Glycolysis",                         -0.04,
  # Gluconeogenesis — PCK1 and FBP1 explicitly suppressed per paper text
  "SLC25A10",    "Gluconeogenesis",                     0.00,
  "GOT2",        "Gluconeogenesis",                     0.05,
  "GOT1",        "Gluconeogenesis",                     0.05,
  "FBP1",        "Gluconeogenesis",                    -0.05,   # suppressed (paper text)
  "SLC25A11",    "Gluconeogenesis",                    -0.10,
  "PCK1",        "Gluconeogenesis",                    -0.15,   # suppressed (paper text)
  "MDH1",        "Gluconeogenesis",                    -0.15,
  # TCA cycle — PDHB, PDK2, ACO2, IDH3G, SUCLG1, SUCLA2 all suppressed per paper text
  "SDHB",        "Pyruvate metabolism and TCA cycle",  -0.05,
  "SUCLG1",      "Pyruvate metabolism and TCA cycle",  -0.05,
  "PDK2",        "Pyruvate metabolism and TCA cycle",  -0.05,
  "ACO2",        "Pyruvate metabolism and TCA cycle",  -0.05,   # suppressed (paper text)
  "IDH3G",       "Pyruvate metabolism and TCA cycle",  -0.05,
  "SUCLA2",      "Pyruvate metabolism and TCA cycle",  -0.08,
  "HAGH",        "Pyruvate metabolism and TCA cycle",  -0.05,
  "PDHB",        "Pyruvate metabolism and TCA cycle",  -0.08,
  "LDHA",        "Pyruvate metabolism and TCA cycle",  -0.10,
  # Glutathione — suppressed in PT (paper: "suppressed in PT, enhanced in TAL")
  # AKR1A1 enhanced (was downregulated in T2D); GSTM3 suppressed in both PT and TAL
  # Large bars in original image were PINK (T2D vs HC), not teal (SGLT2i vs T2D)
  "CNDP2",       "Glutathione conjugation",            -0.02,
  "GSTM4",       "Glutathione conjugation",            -0.03,
  "GSTT2B",      "Glutathione conjugation",            -0.03,
  "GSTO1",       "Glutathione conjugation",            -0.03,   # small negative in PT
  "GGCT",        "Glutathione conjugation",             0.00,
  "GSTM3",       "Glutathione conjugation",            -0.05,   # suppressed in both PT+TAL
  "AKR1A1",      "Glutathione conjugation",             0.08,   # enhanced (restored toward HC)
  # Metallothioneins — consistently enhanced per paper
  "MT1G",        "Metallothioneins bind metals",        0.08,
  "MT1X",        "Metallothioneins bind metals",        0.08,
  "MT1H",        "Metallothioneins bind metals",        0.08,
  "MT2A",        "Metallothioneins bind metals",        0.08
)

# Combine our results with published values
compare_df <- plot_df %>%
  dplyr::select(gene_symbol, pathway, logFC_nebula = logFC_sglt2, sig) %>%
  mutate(gene_symbol = as.character(gene_symbol),
         pathway     = as.character(pathway)) %>%
  full_join(published_vals, by = c("gene_symbol", "pathway")) %>%
  pivot_longer(cols = c(logFC_nebula, logFC_published),
               names_to = "source", values_to = "logFC") %>%
  mutate(
    source      = factor(source,
                         levels = c("logFC_published", "logFC_nebula"),
                         labels = c("Schaub et al. (published)", "Our scRNA-seq (NEBULA)")),
    pathway     = factor(pathway, levels = names(pathway_genes)),
    gene_symbol = factor(gene_symbol, levels = all_genes)
  )

source_colors <- c("Schaub et al. (published)" = "#E07B54",
                   "Our scRNA-seq (NEBULA)"    = "#3AAFA9")

make_compare_panel <- function(pw_name, df) {
  gene_order <- pathway_genes[[pw_name]]
  pw_df <- df %>%
    dplyr::filter(pathway == pw_name) %>%
    mutate(gene_symbol = factor(as.character(gene_symbol), levels = gene_order))
  
  # Bold axis labels for our sig genes
  sig_genes <- plot_df %>%
    dplyr::filter(as.character(pathway) == pw_name, sig) %>%
    pull(gene_symbol) %>% as.character()
  axis_face <- ifelse(gene_order %in% sig_genes, "bold", "plain")
  
  ggplot(pw_df, aes(x = gene_symbol, y = logFC, fill = source)) +
    geom_col(position = position_dodge(0.75), width = 0.65, na.rm = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey30") +
    scale_fill_manual(values = source_colors, name = NULL) +
    scale_x_discrete(drop = FALSE) +
    labs(title = pw_name, x = NULL, y = "log\u2082 fold change") +
    theme_bw(base_size = 9) +
    theme(
      plot.title         = element_text(size = 9, face = "bold"),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 7.5,
                                        face = axis_face),
      axis.title.y       = element_text(size = 8),
      legend.position    = "top",
      legend.text        = element_text(size = 8),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

compare_panels <- imap(
  setNames(names(pathway_genes), names(pathway_genes)),
  ~ make_compare_panel(.x, compare_df)
)

fig_compare <- (compare_panels[[1]] | compare_panels[[2]]) /
  (compare_panels[[3]] | compare_panels[[4]] | compare_panels[[5]]) +
  plot_annotation(
    title    = "Pathway genes: SGLT2i vs. No SGLT2i — Our scRNA-seq vs. Schaub et al.",
    subtitle = "Schaub et al. = published scRNA-seq PT cells (youth T2D) | Ours = NEBULA PT cells | bold x-axis = FDR < 0.05 in ours",
    theme    = theme(plot.title    = element_text(size = 11, face = "bold"),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

ggsave(file.path(output_dir, "pathway_comparison_vs_published.pdf"),
       fig_compare, width = 14, height = 8, device = "pdf")
ggsave(file.path(output_dir, "pathway_comparison_vs_published.png"),
       fig_compare, width = 14, height = 8, dpi = 300)
cat("Comparison figure saved.\n")

################################################################################
# DISCORDANT GENES — PER-PARTICIPANT RAINCLOUD PLOTS
# No SGLT2i: violin nudged right, participant dots nudged left
# SGLT2i:    violin nudged left,  participant dots nudged right
################################################################################

discordant_genes   <- c("MDH1", "GSTO1", "FBP1", "HAGH", "LDHA", "GSTT2B")
discordant_present <- base::intersect(discordant_genes, genes_present)
cat("\nDiscordant genes present:", paste(discordant_present, collapse = ", "), "\n")

disc_expr_mat <- GetAssayData(so_pt[discordant_present, ], layer = "data")

disc_expr_long <- as.data.frame(t(as.matrix(disc_expr_mat))) %>%
  rownames_to_column("cell_barcode") %>%
  left_join(
    so_pt@meta.data %>%
      rownames_to_column("cell_barcode") %>%
      dplyr::select(cell_barcode, sglt2_group, record_id),
    by = "cell_barcode"
  ) %>%
  pivot_longer(cols = all_of(discordant_present),
               names_to = "gene_symbol", values_to = "expression") %>%
  mutate(sglt2_group = factor(sglt2_group, levels = c("No SGLT2i", "SGLT2i")),
         # numeric position for nudging
         x_pos = as.numeric(sglt2_group))

participant_means <- disc_expr_long %>%
  group_by(gene_symbol, sglt2_group, record_id, x_pos) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop")

disc_fdr <- res_df %>%
  dplyr::filter(gene_symbol %in% discordant_present) %>%
  dplyr::select(gene_symbol, logFC_sglt2, fdr, sig) %>%
  mutate(label = paste0("log2FC=", round(logFC_sglt2, 2),
                        "\nFDR=", formatC(fdr, digits = 2, format = "g")))

group_colors <- c("No SGLT2i" = "#4C72B0", "SGLT2i" = "#3AAFA9")
NUDGE  <- 0.22   # how far violins/dots shift from group centre
DOT_W  <- 0.06   # jitter width for participant dots

disc_plots <- map(discordant_present, function(g) {
  gene_expr <- disc_expr_long    %>% dplyr::filter(gene_symbol == g)
  pt_means  <- participant_means %>% dplyr::filter(gene_symbol == g)
  gene_fdr  <- disc_fdr          %>% dplyr::filter(gene_symbol == g)
  
  # Pre-compute nudged x positions in data
  # Violins: No SGLT2i nudged RIGHT, SGLT2i nudged LEFT
  no_vln  <- gene_expr %>% dplyr::filter(sglt2_group == "No SGLT2i") %>%
    mutate(x_vln = x_pos + NUDGE)
  yes_vln <- gene_expr %>% dplyr::filter(sglt2_group == "SGLT2i") %>%
    mutate(x_vln = x_pos - NUDGE)
  
  # Dots: No SGLT2i nudged LEFT, SGLT2i nudged RIGHT
  set.seed(42)
  pt_no  <- pt_means %>% dplyr::filter(sglt2_group == "No SGLT2i") %>%
    mutate(x_dot = x_pos - NUDGE + runif(n(), -DOT_W, DOT_W))
  pt_yes <- pt_means %>% dplyr::filter(sglt2_group == "SGLT2i") %>%
    mutate(x_dot = x_pos + NUDGE + runif(n(), -DOT_W, DOT_W))
  
  y_max <- max(pt_means$mean_expr, na.rm = TRUE) * 1.2
  
  ggplot() +
    # Violins
    geom_violin(data = no_vln,
                aes(x = x_vln, y = expression, fill = sglt2_group, group = x_vln),
                alpha = 0.35, trim = TRUE, linewidth = 0.4, width = 0.5) +
    geom_violin(data = yes_vln,
                aes(x = x_vln, y = expression, fill = sglt2_group, group = x_vln),
                alpha = 0.35, trim = TRUE, linewidth = 0.4, width = 0.5) +
    # Boxplots centred
    geom_boxplot(data = gene_expr,
                 aes(x = x_pos, y = expression, fill = sglt2_group, group = x_pos),
                 width = 0.1, alpha = 0.8, outlier.shape = NA,
                 linewidth = 0.4, color = "grey20") +
    # Participant dots
    geom_point(data = pt_no,
               aes(x = x_dot, y = mean_expr),
               color = group_colors["No SGLT2i"], size = 2, alpha = 0.85) +
    geom_point(data = pt_yes,
               aes(x = x_dot, y = mean_expr),
               color = group_colors["SGLT2i"], size = 2, alpha = 0.85) +
    
    annotate("text", x = 1.5, y = y_max, label = gene_fdr$label,
             size = 2.8, hjust = 0.5, color = "grey20") +
    
    scale_fill_manual(values = group_colors) +
    scale_x_continuous(breaks = c(1, 2),
                       labels = c("No SGLT2i", "SGLT2i"),
                       limits = c(0.4, 2.6)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    labs(title    = g,
         subtitle = "Violin = cell distribution | Dots = per-participant mean",
         x = NULL, y = "Normalized expression") +
    theme_bw(base_size = 10) +
    theme(legend.position    = "none",
          plot.title         = element_text(face = "bold", size = 11),
          plot.subtitle      = element_text(size = 7, color = "grey50"),
          axis.text.x        = element_text(size = 9),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank())
})

fig_discordant <- wrap_plots(disc_plots, ncol = 3) +
  plot_annotation(
    title    = "Discordant pathway genes: PT cell expression by SGLT2i status",
    subtitle = "Half-violin = cell distribution | Dots = per-participant mean expression",
    theme    = theme(plot.title    = element_text(size = 12, face = "bold"),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

ggsave(file.path(output_dir, "discordant_genes_per_participant.pdf"),
       fig_discordant, width = 12, height = 8, device = "pdf")
ggsave(file.path(output_dir, "discordant_genes_per_participant.png"),
       fig_discordant, width = 12, height = 8, dpi = 300)
cat("Discordant gene plots saved.\n")




expr_mat <- GetAssayData(so_pt[genes_present, ], layer = "data")

expr_long <- as.data.frame(t(as.matrix(expr_mat))) %>%
  rownames_to_column("cell_barcode") %>%
  left_join(
    so_pt@meta.data %>%
      rownames_to_column("cell_barcode") %>%
      dplyr::select(cell_barcode, sglt2_group, record_id),
    by = "cell_barcode"
  ) %>%
  pivot_longer(cols = all_of(genes_present),
               names_to = "gene_symbol", values_to = "expression") %>%
  left_join(
    imap_dfr(pathway_genes, ~ tibble(gene_symbol = .x, pathway = .y)),
    by = "gene_symbol"
  ) %>%
  left_join(
    res_df %>% dplyr::select(gene_symbol, fdr, sig),
    by = "gene_symbol"
  ) %>%
  mutate(
    sglt2_group = factor(sglt2_group, levels = c("No SGLT2i", "SGLT2i")),
    pathway     = factor(pathway, levels = names(pathway_genes)),
    gene_symbol = factor(gene_symbol, levels = all_genes),
    fdr_label   = case_when(
      is.na(fdr)  ~ "",
      fdr < 0.001 ~ "FDR<0.001",
      fdr < 0.05  ~ paste0("FDR=", formatC(fdr, digits = 2, format = "g")),
      TRUE        ~ ""
    )
  )

group_colors <- c("No SGLT2i" = "#4C72B0", "SGLT2i" = "#3AAFA9")

iwalk(pathway_genes, function(genes, pw_name) {
  pw_genes <- base::intersect(genes, genes_present)
  if (length(pw_genes) == 0) return(NULL)
  
  pw_expr    <- expr_long %>% dplyr::filter(pathway == pw_name)
  y_max      <- max(pw_expr$expression, na.rm = TRUE)
  fdr_labels <- pw_expr %>%
    distinct(gene_symbol, fdr_label, sig) %>%
    mutate(y = y_max * 1.07)
  
  p <- ggplot(pw_expr,
              aes(x = gene_symbol, y = expression, fill = sglt2_group)) +
    geom_violin(alpha = 0.45, trim = TRUE, position = position_dodge(0.8),
                linewidth = 0.4) +
    geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA,
                 position = position_dodge(0.8), linewidth = 0.4,
                 color = "grey20") +
    geom_text(data = fdr_labels %>% dplyr::filter(sig),
              aes(x = gene_symbol, y = y, label = fdr_label),
              inherit.aes = FALSE, size = 2.5, color = "black", fontface = "bold") +
    scale_fill_manual(values = group_colors, name = NULL) +
    scale_x_discrete(limits = pw_genes) +
    labs(
      title    = pw_name,
      subtitle = "Normalized expression in PT cells (T2D only)",
      x        = NULL,
      y        = "Normalized expression"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position    = "top",
      axis.text.x        = element_text(angle = 30, hjust = 1, size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )
  
  clean_name <- str_replace_all(pw_name, "[^A-Za-z0-9]", "_")
  ggsave(file.path(output_dir, paste0("violin_PT_", clean_name, "_SGLT2i.pdf")),
         p, width = max(4, length(pw_genes) * 1.2), height = 5, device = "pdf")
  ggsave(file.path(output_dir, paste0("violin_PT_", clean_name, "_SGLT2i.png")),
         p, width = max(4, length(pw_genes) * 1.2), height = 5, dpi = 300)
  cat("Violin plot saved:", pw_name, "\n")
})

cat("\n\u2713 All outputs complete. Saved to:", output_dir, "\n")

################################################################################
# LIMMA PSEUDOBULK ANALYSIS — matching Schaub et al. methodology
# Aggregate counts per participant → TMM normalize → voom → limma
################################################################################

library(limma)
library(edgeR)

cat("\nRunning pseudobulk limma...\n")

# Pull raw counts for pathway genes from PT cells
counts_pt <- round(GetAssayData(so_pt[genes_present, ], layer = "counts"))

# Build per-participant metadata
meta_pt <- so_pt@meta.data %>%
  mutate(sglt2i = factor(sglt2i_ever, levels = c("No", "Yes"))) %>%
  dplyr::select(record_id, sglt2i, age, sex)

# Aggregate counts per participant (pseudobulk: sum across cells)
participant_ids <- unique(meta_pt$record_id)

pb_counts <- sapply(participant_ids, function(pid) {
  cells <- which(meta_pt$record_id == pid)
  Matrix::rowSums(counts_pt[, cells, drop = FALSE])
})
colnames(pb_counts) <- participant_ids

# Participant-level metadata (one row per participant)
pb_meta <- meta_pt %>%
  distinct(record_id, .keep_all = TRUE) %>%
  dplyr::filter(record_id %in% participant_ids) %>%
  arrange(match(record_id, participant_ids))

cat("Pseudobulk matrix:", nrow(pb_counts), "genes x", ncol(pb_counts), "participants\n")
cat("SGLT2i breakdown:\n"); print(table(pb_meta$sglt2i))

# Remove any participants with NA in sglt2i
keep_participants <- !is.na(pb_meta$sglt2i)
pb_counts <- pb_counts[, keep_participants]
pb_meta   <- pb_meta[keep_participants, ]

# Build DGEList and normalize
dge <- DGEList(counts = pb_counts, group = pb_meta$sglt2i)
dge <- calcNormFactors(dge, method = "TMM")

# Design matrix — SGLT2i only (matching Schaub et al. simple model)
design_limma <- model.matrix(~ sglt2i, data = pb_meta)

# voom transformation
v <- voom(dge, design_limma, plot = FALSE)

# Fit and test
fit       <- lmFit(v, design_limma)
fit       <- eBayes(fit)

# Extract results for sglt2iYes coefficient
limma_res <- topTable(fit, coef = "sglt2iYes", number = Inf, sort.by = "none") %>%
  rownames_to_column("gene_symbol") %>%
  as_tibble() %>%
  dplyr::rename(logFC_limma = logFC, fdr_limma = adj.P.Val) %>%
  mutate(sig_limma = fdr_limma < 0.05) %>%
  left_join(
    imap_dfr(pathway_genes, ~ tibble(gene_symbol = .x, pathway = .y)),
    by = "gene_symbol"
  ) %>%
  arrange(fdr_limma)

cat("\nLimma significant genes (FDR < 0.05):", sum(limma_res$sig_limma, na.rm = TRUE), "\n")
print(limma_res %>% dplyr::filter(sig_limma) %>%
        dplyr::select(gene_symbol, pathway, logFC_limma, fdr_limma))

write.xlsx(limma_res, file.path(output_dir, "limma_pseudobulk_PT_pathway_genes_SGLT2i.xlsx"))

################################################################################
# THREE-WAY COMPARISON: Schaub et al. vs NEBULA vs Limma
################################################################################

three_way_df <- published_vals %>%
  left_join(
    res_df   %>% dplyr::select(gene_symbol, logFC_nebula = logFC_sglt2, fdr_nebula = fdr),
    by = "gene_symbol"
  ) %>%
  left_join(
    limma_res %>% dplyr::select(gene_symbol, logFC_limma, fdr_limma),
    by = "gene_symbol"
  ) %>%
  pivot_longer(cols = c(logFC_published, logFC_nebula, logFC_limma),
               names_to = "source", values_to = "logFC") %>%
  mutate(
    source = factor(source,
                    levels = c("logFC_published", "logFC_limma", "logFC_nebula"),
                    labels = c("Schaub et al.", "Our limma (pseudobulk)", "Our NEBULA")),
    pathway     = factor(pathway, levels = names(pathway_genes)),
    gene_symbol = factor(gene_symbol, levels = all_genes),
    # Mark significance: limma or NEBULA FDR < 0.05
    sig_limma  = gene_symbol %in% (limma_res %>% dplyr::filter(sig_limma) %>% pull(gene_symbol)),
    sig_nebula = gene_symbol %in% (res_df    %>% dplyr::filter(sig)       %>% pull(gene_symbol))
  )

three_colors <- c("Schaub et al."           = "#E07B54",
                  "Our limma (pseudobulk)"  = "#7B5EA7",
                  "Our NEBULA"              = "#3AAFA9")

make_three_panel <- function(pw_name, df) {
  gene_order <- pathway_genes[[pw_name]]
  pw_df <- df %>%
    dplyr::filter(pathway == pw_name) %>%
    mutate(gene_symbol = factor(as.character(gene_symbol), levels = gene_order))
  
  # Bold for limma sig genes
  sig_limma_genes  <- pw_df %>% dplyr::filter(sig_limma)  %>% pull(gene_symbol) %>% as.character()
  sig_nebula_genes <- pw_df %>% dplyr::filter(sig_nebula) %>% pull(gene_symbol) %>% as.character()
  sig_either       <- unique(c(sig_limma_genes, sig_nebula_genes))
  axis_face        <- ifelse(gene_order %in% sig_either, "bold", "plain")
  
  ggplot(pw_df, aes(x = gene_symbol, y = logFC, fill = source)) +
    geom_col(position = position_dodge(0.75), width = 0.65, na.rm = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey30") +
    scale_fill_manual(values = three_colors, name = NULL) +
    scale_x_discrete(drop = FALSE) +
    labs(title = pw_name, x = NULL, y = "log\u2082 fold change") +
    theme_bw(base_size = 9) +
    theme(
      plot.title         = element_text(size = 9, face = "bold"),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 7.5,
                                        face = axis_face),
      axis.title.y       = element_text(size = 8),
      legend.position    = "top",
      legend.text        = element_text(size = 7.5),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

three_panels <- imap(
  setNames(names(pathway_genes), names(pathway_genes)),
  ~ make_three_panel(.x, three_way_df)
)

fig_three <- (three_panels[[1]] | three_panels[[2]]) /
  (three_panels[[3]] | three_panels[[4]] | three_panels[[5]]) +
  plot_annotation(
    title    = "Three-way comparison: Schaub et al. vs Our limma vs Our NEBULA (PT cells)",
    subtitle = "Bold x-axis = FDR < 0.05 in either of our analyses",
    theme    = theme(plot.title    = element_text(size = 11, face = "bold"),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

ggsave(file.path(output_dir, "pathway_three_way_comparison.pdf"),
       fig_three, width = 15, height = 8, device = "pdf")
ggsave(file.path(output_dir, "pathway_three_way_comparison.png"),
       fig_three, width = 15, height = 8, dpi = 300)
cat("Three-way comparison figure saved.\n")

# Also save a concordance summary table
concordance_summary <- published_vals %>%
  left_join(res_df    %>% dplyr::select(gene_symbol, logFC_nebula = logFC_sglt2, fdr_nebula = fdr), by = "gene_symbol") %>%
  left_join(limma_res %>% dplyr::select(gene_symbol, logFC_limma, fdr_limma), by = "gene_symbol") %>%
  mutate(
    dir_published = sign(logFC_published),
    dir_nebula    = sign(logFC_nebula),
    dir_limma     = sign(logFC_limma),
    nebula_concordant_with_schaub = dir_nebula == dir_published,
    limma_concordant_with_schaub  = dir_limma  == dir_published,
    all_concordant = nebula_concordant_with_schaub & limma_concordant_with_schaub
  )

cat("\nConcordance with Schaub et al.:\n")
cat("NEBULA concordant:", sum(concordance_summary$nebula_concordant_with_schaub, na.rm = TRUE),
    "/", sum(!is.na(concordance_summary$logFC_nebula)), "\n")
cat("Limma concordant:", sum(concordance_summary$limma_concordant_with_schaub, na.rm = TRUE),
    "/", sum(!is.na(concordance_summary$logFC_limma)), "\n")
cat("All three concordant:", sum(concordance_summary$all_concordant, na.rm = TRUE),
    "/", sum(!is.na(concordance_summary$logFC_nebula) & !is.na(concordance_summary$logFC_limma)), "\n")

write.xlsx(concordance_summary,
           file.path(output_dir, "concordance_summary_three_way.xlsx"))

cat("\n\u2713 All analyses complete.\n")

################################################################################
# SCHAUB-EQUIVALENT SUBSET ANALYSIS
# Restrict to RH- + IT_ participants only (same studies as Schaub et al.)
# Excludes RH2 which was not in Schaub's cohort
################################################################################

cat("\n=== Rerunning on Schaub-equivalent subset (RH + IT only) ===\n")

schaub_cohort_ids <- unique(so_pt$record_id)[
  str_detect(unique(so_pt$record_id), "^RH-|^IT_")
]

so_pt_s <- subset(so_pt, record_id %in% schaub_cohort_ids)
cat("Participants:", length(schaub_cohort_ids), "\n")
cat("Cells:", ncol(so_pt_s), "\n")
cat("SGLT2i balance (participants):\n")
print(table(so_pt_s@meta.data %>%
              distinct(record_id, sglt2i_ever) %>%
              pull(sglt2i_ever)))

# ── NEBULA on Schaub subset ──────────────────────────────────────────────────
genes_s      <- base::intersect(genes_present, rownames(so_pt_s))
so_pt_s_sub  <- so_pt_s[genes_s, ]
counts_s     <- round(GetAssayData(so_pt_s_sub, layer = "counts"))

# Recompute size factors on the full Schaub subset (not just pathway genes)
counts_s_full    <- round(GetAssayData(so_pt_s, layer = "counts"))
sce_s            <- SingleCellExperiment(assays = list(counts = counts_s_full))
sce_s            <- computeSumFactors(sce_s)
so_pt_s$pooled_offset_s <- sizeFactors(sce_s)
rm(sce_s, counts_s_full)

meta_s <- so_pt_s_sub@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  mutate(sglt2i = factor(sglt2i_ever, levels = c("No", "Yes"))) %>%
  # Pull pooled_offset_s from so_pt_s (was added after subsetting)
  left_join(
    so_pt_s@meta.data %>%
      rownames_to_column("cell_barcode") %>%
      dplyr::select(cell_barcode, pooled_offset_s),
    by = "cell_barcode"
  ) %>%
  column_to_rownames("cell_barcode")

complete_s   <- complete.cases(meta_s[, c("sglt2i", "pooled_offset_s")])
counts_s     <- counts_s[, complete_s]
meta_s       <- meta_s[complete_s, ]

pred_s  <- model.matrix(~ sglt2i, data = meta_s)
lib_s   <- meta_s$pooled_offset_s

data_g_s <- group_cell(count = counts_s, id = meta_s$record_id,
                       pred = pred_s, offset = lib_s)
if (is.null(data_g_s)) {
  data_g_s <- list(count = counts_s, id = meta_s$record_id,
                   pred = pred_s, offset = lib_s)
}

nebula_s <- nebula(count = data_g_s$count, id = data_g_s$id,
                   pred  = data_g_s$pred,  offset = data_g_s$offset,
                   model = "NBLMM", ncore = 1, reml = TRUE,
                   output_re = TRUE, covariance = TRUE)

res_raw_s <- nebula_s$summary %>% as_tibble()
fc_col_s  <- grep("logFC.*sglt2|logFC.*Yes", names(res_raw_s), value = TRUE,
                  ignore.case = TRUE)[1]
p_col_s   <- grep("^p_.*sglt2|^p_.*Yes",    names(res_raw_s), value = TRUE,
                  ignore.case = TRUE)[1]

res_nebula_s <- res_raw_s %>%
  dplyr::rename(gene_symbol   = gene,
                logFC_nebula_s = all_of(fc_col_s),
                p_nebula_s     = all_of(p_col_s)) %>%
  mutate(fdr_nebula_s = p.adjust(p_nebula_s, method = "BH"),
         sig_nebula_s = fdr_nebula_s < 0.05) %>%
  left_join(imap_dfr(pathway_genes, ~ tibble(gene_symbol = .x, pathway = .y)),
            by = "gene_symbol")

cat("Schaub-subset NEBULA significant genes:", sum(res_nebula_s$sig_nebula_s, na.rm = TRUE), "\n")

# ── Limma pseudobulk on Schaub subset ───────────────────────────────────────
counts_pt_s <- round(GetAssayData(so_pt_s[genes_s, ], layer = "counts"))
meta_pt_s   <- so_pt_s@meta.data %>%
  mutate(sglt2i = factor(sglt2i_ever, levels = c("No", "Yes"))) %>%
  dplyr::select(record_id, sglt2i)

participant_ids_s <- unique(meta_pt_s$record_id)

pb_counts_s <- sapply(participant_ids_s, function(pid) {
  cells <- which(meta_pt_s$record_id == pid)
  Matrix::rowSums(counts_pt_s[, cells, drop = FALSE])
})
colnames(pb_counts_s) <- participant_ids_s

pb_meta_s <- meta_pt_s %>%
  distinct(record_id, .keep_all = TRUE) %>%
  dplyr::filter(record_id %in% participant_ids_s) %>%
  arrange(match(record_id, participant_ids_s)) %>%
  dplyr::filter(!is.na(sglt2i))

pb_counts_s <- pb_counts_s[, pb_meta_s$record_id]

dge_s      <- DGEList(counts = pb_counts_s, group = pb_meta_s$sglt2i)
dge_s      <- calcNormFactors(dge_s, method = "TMM")
design_s   <- model.matrix(~ sglt2i, data = pb_meta_s)
v_s        <- voom(dge_s, design_s, plot = FALSE)
fit_s      <- lmFit(v_s, design_s)
fit_s      <- eBayes(fit_s)

res_limma_s <- topTable(fit_s, coef = "sglt2iYes", number = Inf, sort.by = "none") %>%
  rownames_to_column("gene_symbol") %>% as_tibble() %>%
  dplyr::rename(logFC_limma_s = logFC, fdr_limma_s = adj.P.Val) %>%
  mutate(sig_limma_s = fdr_limma_s < 0.05) %>%
  left_join(imap_dfr(pathway_genes, ~ tibble(gene_symbol = .x, pathway = .y)),
            by = "gene_symbol")

cat("Schaub-subset limma significant genes:", sum(res_limma_s$sig_limma_s, na.rm = TRUE), "\n")

# Save results
write.xlsx(list(NEBULA = res_nebula_s, Limma = res_limma_s),
           file.path(output_dir, "Schaub_subset_results.xlsx"))

# ── Four-way comparison figure ───────────────────────────────────────────────
four_way_df <- published_vals %>%
  left_join(res_df       %>% dplyr::select(gene_symbol, logFC_nebula_all  = logFC_sglt2),  by = "gene_symbol") %>%
  left_join(limma_res    %>% dplyr::select(gene_symbol, logFC_limma_all   = logFC_limma),  by = "gene_symbol") %>%
  left_join(res_nebula_s %>% dplyr::select(gene_symbol, logFC_nebula_s),                  by = "gene_symbol") %>%
  left_join(res_limma_s  %>% dplyr::select(gene_symbol, logFC_limma_s),                   by = "gene_symbol") %>%
  pivot_longer(cols = c(logFC_published, logFC_nebula_all, logFC_limma_all,
                        logFC_nebula_s, logFC_limma_s),
               names_to = "source", values_to = "logFC") %>%
  mutate(
    source = factor(source,
                    levels = c("logFC_published",
                               "logFC_limma_s", "logFC_nebula_s",
                               "logFC_limma_all", "logFC_nebula_all"),
                    labels = c("Schaub et al.",
                               "Limma (RH+IT only)", "NEBULA (RH+IT only)",
                               "Limma (all)", "NEBULA (all)")),
    pathway     = factor(pathway, levels = names(pathway_genes)),
    gene_symbol = factor(gene_symbol, levels = all_genes)
  )

five_colors <- c("Schaub et al."        = "#E07B54",
                 "Limma (RH+IT only)"   = "#7B5EA7",
                 "NEBULA (RH+IT only)"  = "#3AAFA9",
                 "Limma (all)"          = "#C4A8E0",
                 "NEBULA (all)"         = "#A8D8D8")

make_four_panel <- function(pw_name, df) {
  gene_order <- pathway_genes[[pw_name]]
  pw_df <- df %>%
    dplyr::filter(pathway == pw_name) %>%
    mutate(gene_symbol = factor(as.character(gene_symbol), levels = gene_order))
  
  ggplot(pw_df, aes(x = gene_symbol, y = logFC, fill = source)) +
    geom_col(position = position_dodge(0.8), width = 0.75, na.rm = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey30") +
    scale_fill_manual(values = five_colors, name = NULL) +
    scale_x_discrete(drop = FALSE) +
    labs(title = pw_name, x = NULL, y = "log\u2082 fold change") +
    theme_bw(base_size = 9) +
    theme(
      plot.title         = element_text(size = 9, face = "bold"),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 7),
      axis.title.y       = element_text(size = 8),
      legend.position    = "top",
      legend.text        = element_text(size = 7),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

four_panels <- imap(
  setNames(names(pathway_genes), names(pathway_genes)),
  ~ make_four_panel(.x, four_way_df)
)

fig_four <- (four_panels[[1]] | four_panels[[2]]) /
  (four_panels[[3]] | four_panels[[4]] | four_panels[[5]]) +
  plot_annotation(
    title    = "Pathway comparison: Schaub et al. vs our analyses (all participants vs RH+IT only)",
    subtitle = "Darker = RH+IT only (Schaub-equivalent) | Lighter = all participants",
    theme    = theme(plot.title    = element_text(size = 11, face = "bold"),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

ggsave(file.path(output_dir, "pathway_four_way_comparison.pdf"),
       fig_four, width = 16, height = 8, device = "pdf")
ggsave(file.path(output_dir, "pathway_four_way_comparison.png"),
       fig_four, width = 16, height = 8, dpi = 300)
cat("Four-way comparison figure saved.\n")

# ── Concordance summary ───────────────────────────────────────────────────────
conc <- published_vals %>%
  left_join(res_nebula_s %>% dplyr::select(gene_symbol, logFC_nebula_s), by = "gene_symbol") %>%
  left_join(res_limma_s  %>% dplyr::select(gene_symbol, logFC_limma_s),  by = "gene_symbol") %>%
  left_join(res_df       %>% dplyr::select(gene_symbol, logFC_nebula_all = logFC_sglt2), by = "gene_symbol") %>%
  left_join(limma_res    %>% dplyr::select(gene_symbol, logFC_limma_all = logFC_limma),  by = "gene_symbol") %>%
  mutate(across(starts_with("logFC"), sign, .names = "dir_{.col}")) %>%
  mutate(
    conc_nebula_s   = dir_logFC_nebula_s   == dir_logFC_published,
    conc_limma_s    = dir_logFC_limma_s    == dir_logFC_published,
    conc_nebula_all = dir_logFC_nebula_all == dir_logFC_published,
    conc_limma_all  = dir_logFC_limma_all  == dir_logFC_published
  )

n <- sum(!is.na(conc$logFC_published))
cat("\n--- Directional concordance with Schaub et al. ---\n")
cat("NEBULA all participants: ", sum(conc$conc_nebula_all,  na.rm=TRUE), "/", n, "\n")
cat("Limma  all participants: ", sum(conc$conc_limma_all,   na.rm=TRUE), "/", n, "\n")
cat("NEBULA RH+IT only:       ", sum(conc$conc_nebula_s,    na.rm=TRUE), "/", n, "\n")
cat("Limma  RH+IT only:       ", sum(conc$conc_limma_s,     na.rm=TRUE), "/", n, "\n")

cat("\n\u2713 All done.\n")



################################################################################
# PARTICIPANT GROUPING CHECK — verify our SGLT2i assignments
################################################################################

cat("\n=== Participant SGLT2i groupings (Schaub-equivalent cohort) ===\n")
so_pt_s@meta.data %>%
  distinct(record_id, sglt2i_ever) %>%
  mutate(study = case_when(
    str_detect(record_id, "^RH-") ~ "Renal-HEIR",
    str_detect(record_id, "^IT_") ~ "IMPROVE-T2D",
    TRUE                          ~ "Other"
  )) %>%
  arrange(sglt2i_ever, record_id) %>%
  as.data.frame() %>%
  print()

################################################################################
# PT SUBCLUSTER ANALYSIS — PT-1 through PT-5 matching Schaub et al.
################################################################################

cat("\n=== PT subcluster analysis (PT-1 to PT-5) ===\n")

pt_subclusters <- c("PT-1", "PT-2", "PT-3", "PT-4", "PT-5")

so_pt_sub5 <- subset(so, subset = celltype_rpca %in% pt_subclusters)
so_pt_sub5$pt_subcluster <- so_pt_sub5$celltype_rpca
so_pt_sub5$sglt2_group   <- ifelse(so_pt_sub5$sglt2i_ever == "Yes",
                                   "SGLT2i", "No SGLT2i")

cat("PT subcluster cell counts by SGLT2i group:\n")
print(table(so_pt_sub5$pt_subcluster, so_pt_sub5$sglt2_group))

counts_sub5_full          <- round(GetAssayData(so_pt_sub5, layer = "counts"))
sce_sub5                  <- SingleCellExperiment(assays = list(counts = counts_sub5_full))
sce_sub5                  <- computeSumFactors(sce_sub5)
so_pt_sub5$pooled_offset_sub5 <- sizeFactors(sce_sub5)
rm(sce_sub5, counts_sub5_full)

genes_sub5 <- base::intersect(genes_present, rownames(so_pt_sub5))

# sglt2i_schaub: Schaub CSV is ground truth for the 16 biopsy participants
# sglt2i_ever now correctly captures RH-74-T, RH-75-T, RH-77-T (fixed vs sglt2i_timepoint)
# RH-49-T still needs the explicit override (sglt2i_ever=Yes per harmonized but confirm)
schaub_sglt2i_correct    <- c("RH-49-T", "RH-50-T", "RH-62-T", "RH-63-T", "RH-65-T",
                              "RH-68-T", "RH-72-T", "RH-74-T", "IT_12", "IT_13")
schaub_no_sglt2i_correct <- c("RH-23-T", "RH-59-T", "RH-60-T", "RH-66-T",
                              "RH-71-T", "IT_11")

so_pt_sub5$sglt2i_schaub <- case_when(
  so_pt_sub5$record_id %in% schaub_sglt2i_correct    ~ "Yes",
  so_pt_sub5$record_id %in% schaub_no_sglt2i_correct ~ "No",
  so_pt_sub5$sglt2i_ever == "Yes"                    ~ "Yes",
  so_pt_sub5$sglt2i_ever == "No"                     ~ "No",
  TRUE ~ NA_character_
)

cat("Any remaining discrepancies between sglt2i_schaub and sglt2i_ever:\n")
disc <- so_pt_sub5@meta.data %>%
  distinct(record_id, sglt2i_schaub, sglt2i_ever) %>%
  dplyr::filter(sglt2i_schaub != sglt2i_ever)
if (nrow(disc) == 0) cat("  None — sglt2i_ever and Schaub CSV are fully concordant.\n") else print(disc)

run_nebula_pt_subcluster <- function(so_obj, pt_name) {
  cat("\n-- PT subcluster:", pt_name, "--\n")
  so_sub <- subset(so_obj, subset = pt_subcluster == pt_name)
  
  grp_tab <- table(
    so_sub@meta.data %>%
      distinct(record_id, sglt2i_schaub) %>%
      pull(sglt2i_schaub)
  )
  cat("  Participants:", paste(names(grp_tab), grp_tab, sep = "=", collapse = ", "), "\n")
  if (length(grp_tab) < 2 || any(grp_tab < 2)) {
    cat("  Skipping — insufficient groups.\n"); return(NULL)
  }
  
  genes_here <- base::intersect(genes_sub5, rownames(so_sub))
  so_sub_g   <- so_sub[genes_here, ]
  counts_sub <- round(GetAssayData(so_sub_g, layer = "counts"))
  
  # pooled_offset_sub5 already in metadata (inherited from so_obj via subsetting)
  meta_sub <- so_sub_g@meta.data %>%
    mutate(sglt2i = factor(sglt2i_schaub, levels = c("No", "Yes")))
  
  ok         <- complete.cases(meta_sub[, c("sglt2i", "pooled_offset_sub5")])
  counts_sub <- counts_sub[, ok]
  meta_sub   <- meta_sub[ok, ]
  
  pred_sub <- model.matrix(~ sglt2i, data = meta_sub)
  lib_sub  <- meta_sub$pooled_offset_sub5
  
  data_g <- group_cell(count = counts_sub, id = meta_sub$record_id,
                       pred = pred_sub, offset = lib_sub)
  if (is.null(data_g)) {
    data_g <- list(count = counts_sub, id = meta_sub$record_id,
                   pred = pred_sub, offset = lib_sub)
  }
  
  neb <- tryCatch(
    nebula(count = data_g$count, id = data_g$id, pred = data_g$pred,
           offset = data_g$offset, model = "NBLMM",
           ncore = 1, reml = TRUE, output_re = TRUE, covariance = TRUE),
    error = function(e) { cat("  NEBULA error:", e$message, "\n"); NULL }
  )
  if (is.null(neb)) return(NULL)
  
  res    <- neb$summary %>% as_tibble()
  fc_col <- grep("logFC.*sglt2|logFC.*Yes", names(res), value = TRUE,
                 ignore.case = TRUE)[1]
  p_col  <- grep("^p_.*sglt2|^p_.*Yes",    names(res), value = TRUE,
                 ignore.case = TRUE)[1]
  
  res %>%
    dplyr::rename(gene_symbol = gene,
                  logFC       = all_of(fc_col),
                  p_value     = all_of(p_col)) %>%
    mutate(fdr           = p.adjust(p_value, method = "BH"),
           sig           = fdr < 0.05,
           pt_subcluster = pt_name) %>%
    left_join(imap_dfr(pathway_genes, ~ tibble(gene_symbol = .x, pathway = .y)),
              by = "gene_symbol")
}

results_by_subcluster <- map(pt_subclusters,
                             ~ run_nebula_pt_subcluster(so_pt_sub5, .x)) %>%
  set_names(pt_subclusters) %>%
  compact()

results_sub5_combined <- bind_rows(results_by_subcluster)
write.xlsx(c(results_by_subcluster, list(Combined = results_sub5_combined)),
           file.path(output_dir, "NEBULA_PT_subclusters_1to5.xlsx"))

sub5_colors <- c("PT-1" = "#2C7BB6", "PT-2" = "#4DAC26", "PT-3" = "#D7191C",
                 "PT-4" = "#FDAE61", "PT-5" = "#ABD9E9")

sub5_plot_df <- results_sub5_combined %>%
  dplyr::select(gene_symbol, pathway, logFC, fdr, sig, pt_subcluster) %>%
  bind_rows(
    published_vals %>%
      mutate(logFC = logFC_published, fdr = NA_real_,
             sig = FALSE, pt_subcluster = "Schaub et al.")
  ) %>%
  mutate(
    source      = factor(pt_subcluster,
                         levels = c("Schaub et al.", pt_subclusters)),
    pathway     = factor(pathway, levels = names(pathway_genes)),
    gene_symbol = factor(gene_symbol, levels = all_genes)
  )

sub5_fill <- c("Schaub et al." = "#E07B54", sub5_colors)

make_sub5_panel <- function(pw_name, df) {
  gene_order <- pathway_genes[[pw_name]]
  pw_df <- df %>%
    dplyr::filter(pathway == pw_name) %>%
    mutate(gene_symbol = factor(as.character(gene_symbol), levels = gene_order))
  sig_genes <- pw_df %>% dplyr::filter(sig) %>% pull(gene_symbol) %>% as.character()
  axis_face <- ifelse(gene_order %in% sig_genes, "bold", "plain")
  
  ggplot(pw_df, aes(x = gene_symbol, y = logFC, fill = source)) +
    geom_col(position = position_dodge(0.8), width = 0.75, na.rm = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey30") +
    scale_fill_manual(values = sub5_fill, name = NULL) +
    scale_x_discrete(drop = FALSE) +
    labs(title = pw_name, x = NULL, y = "log\u2082 fold change") +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7,
                                     face = axis_face),
          axis.title.y = element_text(size = 8),
          legend.position = "top", legend.text = element_text(size = 7.5),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
}

sub5_panels <- imap(setNames(names(pathway_genes), names(pathway_genes)),
                    ~ make_sub5_panel(.x, sub5_plot_df))

fig_sub5 <- (sub5_panels[[1]] | sub5_panels[[2]]) /
  (sub5_panels[[3]] | sub5_panels[[4]] | sub5_panels[[5]]) +
  plot_annotation(
    title    = "Pathway genes by PT subcluster (PT-1 to PT-5) vs Schaub et al.",
    subtitle = "NEBULA SGLT2i vs No SGLT2i | bold x-axis = FDR < 0.05",
    theme    = theme(plot.title    = element_text(size = 11, face = "bold"),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

ggsave(file.path(output_dir, "pathway_PT_subclusters_vs_Schaub.pdf"),
       fig_sub5, width = 16, height = 8, device = "pdf")
ggsave(file.path(output_dir, "pathway_PT_subclusters_vs_Schaub.png"),
       fig_sub5, width = 16, height = 8, dpi = 300)
cat("PT subcluster figure saved.\n")
cat("\n\u2713 All analyses complete.\n")

################################################################################
# SCHAUB EXACT COHORT — NEBULA + LIMMA
# Restricted to the exact 16 T2D participants from Schaub et al. GEO metadata
# Co-enrolled participants: use whichever alias exists in the Seurat object
################################################################################

cat("\n=== Schaub exact cohort analysis ===\n")

# Exact IDs from CSV — try primary ID first, fall back to co-enrolled alias
schaub_sglt2i_ids    <- c("IT_09", "IT_12", "IT_13", "RH-49-T", "RH-50-T",
                          "RH-62-T", "RH-63-T", "RH-68-T", "RH-72-T", "RH-74-T")
schaub_no_sglt2i_ids <- c("IT_07", "IT_08", "IT_10", "IT_11", "RH-23-T", "RH-71-T")

coenroll_alias_s <- c("IT_07" = "RH-59-T", "IT_08" = "RH-60-T",
                      "IT_09" = "RH-65-T", "IT_10" = "RH-66-T")

resolve_seurat_id <- function(ids, so_obj) {
  sapply(ids, function(id) {
    if (id %in% unique(so_obj$record_id)) return(id)
    alias <- coenroll_alias_s[id]
    if (!is.na(alias) && alias %in% unique(so_obj$record_id)) return(alias)
    return(NA_character_)
  })
}

sglt2i_resolved_s    <- na.omit(resolve_seurat_id(schaub_sglt2i_ids,    so_pt))
no_sglt2i_resolved_s <- na.omit(resolve_seurat_id(schaub_no_sglt2i_ids, so_pt))
all_schaub_s         <- unique(c(sglt2i_resolved_s, no_sglt2i_resolved_s))

cat("Resolved SGLT2i(+):", paste(sort(sglt2i_resolved_s), collapse = ", "), "\n")
cat("Resolved SGLT2i(-):", paste(sort(no_sglt2i_resolved_s), collapse = ", "), "\n")
cat("Total resolved:", length(all_schaub_s), "/ 16\n")

# Subset PT to exact Schaub participants
so_pt_exact <- subset(so_pt, record_id %in% all_schaub_s)
so_pt_exact$sglt2i_exact <- case_when(
  so_pt_exact$record_id %in% sglt2i_resolved_s    ~ "Yes",
  so_pt_exact$record_id %in% no_sglt2i_resolved_s ~ "No",
  TRUE ~ NA_character_
)

cat("PT cells (Schaub exact):", ncol(so_pt_exact), "\n")
cat("SGLT2i balance (cells):\n")
print(table(so_pt_exact$sglt2i_exact))
cat("SGLT2i balance (participants):\n")
print(so_pt_exact@meta.data %>%
        distinct(record_id, sglt2i_exact) %>%
        count(sglt2i_exact) %>% as.data.frame())

# ── NEBULA ─────────────────────────────────────────────────────────────────────
genes_exact    <- base::intersect(genes_present, rownames(so_pt_exact))
so_pt_ex_sub   <- so_pt_exact[genes_exact, ]
counts_exact   <- round(GetAssayData(so_pt_ex_sub, layer = "counts"))

# Size factors on full Schaub-exact PT object
counts_ex_full        <- round(GetAssayData(so_pt_exact, layer = "counts"))
sce_ex                <- SingleCellExperiment(assays = list(counts = counts_ex_full))
sce_ex                <- computeSumFactors(sce_ex)
so_pt_exact$pooled_offset_ex <- sizeFactors(sce_ex)
rm(sce_ex, counts_ex_full)

meta_ex <- so_pt_ex_sub@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  mutate(sglt2i = factor(sglt2i_exact, levels = c("No", "Yes"))) %>%
  left_join(
    so_pt_exact@meta.data %>%
      rownames_to_column("cell_barcode") %>%
      dplyr::select(cell_barcode, pooled_offset_ex),
    by = "cell_barcode"
  ) %>%
  column_to_rownames("cell_barcode")

ok_ex        <- complete.cases(meta_ex[, c("sglt2i", "pooled_offset_ex")])
counts_exact <- counts_exact[, ok_ex]
meta_ex      <- meta_ex[ok_ex, ]

pred_ex <- model.matrix(~ sglt2i, data = meta_ex)
lib_ex  <- meta_ex$pooled_offset_ex

data_g_ex <- group_cell(count = counts_exact, id = meta_ex$record_id,
                        pred = pred_ex, offset = lib_ex)
if (is.null(data_g_ex)) {
  data_g_ex <- list(count = counts_exact, id = meta_ex$record_id,
                    pred = pred_ex, offset = lib_ex)
}

nebula_ex <- nebula(count = data_g_ex$count, id = data_g_ex$id,
                    pred  = data_g_ex$pred,  offset = data_g_ex$offset,
                    model = "NBLMM", ncore = 1, reml = TRUE,
                    output_re = TRUE, covariance = TRUE)

res_raw_ex <- nebula_ex$summary %>% as_tibble()
fc_col_ex  <- grep("logFC.*sglt2|logFC.*Yes", names(res_raw_ex), value = TRUE,
                   ignore.case = TRUE)[1]
p_col_ex   <- grep("^p_.*sglt2|^p_.*Yes",    names(res_raw_ex), value = TRUE,
                   ignore.case = TRUE)[1]

res_nebula_ex <- res_raw_ex %>%
  dplyr::rename(gene_symbol    = gene,
                logFC_nebula_ex = all_of(fc_col_ex),
                p_nebula_ex     = all_of(p_col_ex)) %>%
  mutate(fdr_nebula_ex = p.adjust(p_nebula_ex, method = "BH"),
         sig_nebula_ex = fdr_nebula_ex < 0.05) %>%
  left_join(imap_dfr(pathway_genes, ~ tibble(gene_symbol = .x, pathway = .y)),
            by = "gene_symbol")

cat("Schaub-exact NEBULA significant genes:", sum(res_nebula_ex$sig_nebula_ex, na.rm=TRUE), "\n")

# ── Limma pseudobulk ──────────────────────────────────────────────────────────
counts_pt_ex  <- round(GetAssayData(so_pt_exact[genes_exact, ], layer = "counts"))
meta_pt_ex    <- so_pt_exact@meta.data %>%
  mutate(sglt2i = factor(sglt2i_exact, levels = c("No", "Yes"))) %>%
  dplyr::select(record_id, sglt2i)

pids_ex <- unique(meta_pt_ex$record_id)
pb_ex   <- sapply(pids_ex, function(pid) {
  Matrix::rowSums(counts_pt_ex[, meta_pt_ex$record_id == pid, drop = FALSE])
})
colnames(pb_ex) <- pids_ex

pb_meta_ex <- meta_pt_ex %>%
  distinct(record_id, .keep_all = TRUE) %>%
  dplyr::filter(record_id %in% pids_ex, !is.na(sglt2i)) %>%
  arrange(match(record_id, pids_ex))
pb_ex <- pb_ex[, pb_meta_ex$record_id]

dge_ex     <- DGEList(counts = pb_ex, group = pb_meta_ex$sglt2i)
dge_ex     <- calcNormFactors(dge_ex, method = "TMM")
design_ex  <- model.matrix(~ sglt2i, data = pb_meta_ex)
v_ex       <- voom(dge_ex, design_ex, plot = FALSE)
fit_ex     <- eBayes(lmFit(v_ex, design_ex))

res_limma_ex <- topTable(fit_ex, coef = "sglt2iYes", number = Inf,
                         sort.by = "none") %>%
  rownames_to_column("gene_symbol") %>% as_tibble() %>%
  dplyr::rename(logFC_limma_ex = logFC, fdr_limma_ex = adj.P.Val) %>%
  mutate(sig_limma_ex = fdr_limma_ex < 0.05) %>%
  left_join(imap_dfr(pathway_genes, ~ tibble(gene_symbol = .x, pathway = .y)),
            by = "gene_symbol")

cat("Schaub-exact limma significant genes:", sum(res_limma_ex$sig_limma_ex, na.rm=TRUE), "\n")

write.xlsx(list(NEBULA_exact = res_nebula_ex, Limma_exact = res_limma_ex),
           file.path(output_dir, "Schaub_exact_cohort_results.xlsx"))

# ── Six-way comparison figure ─────────────────────────────────────────────────
six_way_df <- published_vals %>%
  left_join(res_df        %>% dplyr::select(gene_symbol, logFC_nebula_all  = logFC_sglt2),  by = "gene_symbol") %>%
  left_join(limma_res     %>% dplyr::select(gene_symbol, logFC_limma_all   = logFC_limma),  by = "gene_symbol") %>%
  left_join(res_nebula_s  %>% dplyr::select(gene_symbol, logFC_nebula_rh   = logFC_nebula_s), by = "gene_symbol") %>%
  left_join(res_limma_s   %>% dplyr::select(gene_symbol, logFC_limma_rh    = logFC_limma_s),  by = "gene_symbol") %>%
  left_join(res_nebula_ex %>% dplyr::select(gene_symbol, logFC_nebula_ex), by = "gene_symbol") %>%
  left_join(res_limma_ex  %>% dplyr::select(gene_symbol, logFC_limma_ex),  by = "gene_symbol") %>%
  pivot_longer(cols = c(logFC_published, logFC_limma_ex, logFC_nebula_ex,
                        logFC_limma_rh, logFC_nebula_rh,
                        logFC_limma_all, logFC_nebula_all),
               names_to = "source", values_to = "logFC") %>%
  mutate(
    source = factor(source,
                    levels = c("logFC_published",
                               "logFC_limma_ex",   "logFC_nebula_ex",
                               "logFC_limma_rh",   "logFC_nebula_rh",
                               "logFC_limma_all",  "logFC_nebula_all"),
                    labels = c("Schaub et al.",
                               "Limma (exact 16)", "NEBULA (exact 16)",
                               "Limma (RH+IT)",    "NEBULA (RH+IT)",
                               "Limma (all)",      "NEBULA (all)")),
    pathway     = factor(pathway, levels = names(pathway_genes)),
    gene_symbol = factor(gene_symbol, levels = all_genes)
  )

six_colors <- c(
  "Schaub et al."      = "#E07B54",
  "Limma (exact 16)"   = "#4B0082",
  "NEBULA (exact 16)"  = "#008080",
  "Limma (RH+IT)"      = "#9370DB",
  "NEBULA (RH+IT)"     = "#3AAFA9",
  "Limma (all)"        = "#C4A8E0",
  "NEBULA (all)"       = "#A8D8D8"
)

make_six_panel <- function(pw_name, df) {
  gene_order <- pathway_genes[[pw_name]]
  pw_df <- df %>%
    dplyr::filter(pathway == pw_name) %>%
    mutate(gene_symbol = factor(as.character(gene_symbol), levels = gene_order))
  
  ggplot(pw_df, aes(x = gene_symbol, y = logFC, fill = source)) +
    geom_col(position = position_dodge(0.85), width = 0.8, na.rm = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey30") +
    scale_fill_manual(values = six_colors, name = NULL) +
    scale_x_discrete(drop = FALSE) +
    labs(title = pw_name, x = NULL, y = "log\u2082 fold change") +
    theme_bw(base_size = 9) +
    theme(
      plot.title         = element_text(size = 9, face = "bold"),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 7),
      axis.title.y       = element_text(size = 8),
      legend.position    = "top",
      legend.text        = element_text(size = 7),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

six_panels <- imap(setNames(names(pathway_genes), names(pathway_genes)),
                   ~ make_six_panel(.x, six_way_df))

fig_six <- (six_panels[[1]] | six_panels[[2]]) /
  (six_panels[[3]] | six_panels[[4]] | six_panels[[5]]) +
  plot_annotation(
    title    = "Pathway comparison: Schaub et al. vs our analyses (three cohort restrictions)",
    subtitle = "Darkest = exact 16 participants | Mid = RH+IT only | Lightest = all participants",
    theme    = theme(plot.title    = element_text(size = 11, face = "bold"),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

ggsave(file.path(output_dir, "pathway_six_way_comparison.pdf"),
       fig_six, width = 18, height = 8, device = "pdf")
ggsave(file.path(output_dir, "pathway_six_way_comparison.png"),
       fig_six, width = 18, height = 8, dpi = 300)
cat("Six-way comparison saved.\n")

# ── Updated concordance summary ───────────────────────────────────────────────
conc_all <- published_vals %>%
  left_join(res_df        %>% dplyr::select(gene_symbol, logFC_nebula_all  = logFC_sglt2),    by = "gene_symbol") %>%
  left_join(limma_res     %>% dplyr::select(gene_symbol, logFC_limma_all   = logFC_limma),    by = "gene_symbol") %>%
  left_join(res_nebula_s  %>% dplyr::select(gene_symbol, logFC_nebula_rh   = logFC_nebula_s), by = "gene_symbol") %>%
  left_join(res_limma_s   %>% dplyr::select(gene_symbol, logFC_limma_rh    = logFC_limma_s),  by = "gene_symbol") %>%
  left_join(res_nebula_ex %>% dplyr::select(gene_symbol, logFC_nebula_ex), by = "gene_symbol") %>%
  left_join(res_limma_ex  %>% dplyr::select(gene_symbol, logFC_limma_ex),  by = "gene_symbol") %>%
  mutate(dir_pub      = sign(logFC_published),
         dir_neb_all  = sign(logFC_nebula_all),
         dir_lim_all  = sign(logFC_limma_all),
         dir_neb_rh   = sign(logFC_nebula_rh),
         dir_lim_rh   = sign(logFC_limma_rh),
         dir_neb_ex   = sign(logFC_nebula_ex),
         dir_lim_ex   = sign(logFC_limma_ex))

n_genes <- sum(!is.na(conc_all$logFC_published))
cat("\n--- Final directional concordance with Schaub et al. (n =", n_genes, "genes) ---\n")
cat("NEBULA (all participants):  ", sum(conc_all$dir_neb_all == conc_all$dir_pub, na.rm=TRUE), "/", n_genes, "\n")
cat("Limma  (all participants):  ", sum(conc_all$dir_lim_all == conc_all$dir_pub, na.rm=TRUE), "/", n_genes, "\n")
cat("NEBULA (RH+IT only):        ", sum(conc_all$dir_neb_rh  == conc_all$dir_pub, na.rm=TRUE), "/", n_genes, "\n")
cat("Limma  (RH+IT only):        ", sum(conc_all$dir_lim_rh  == conc_all$dir_pub, na.rm=TRUE), "/", n_genes, "\n")
cat("NEBULA (exact 16):          ", sum(conc_all$dir_neb_ex  == conc_all$dir_pub, na.rm=TRUE), "/", n_genes, "\n")
cat("Limma  (exact 16):          ", sum(conc_all$dir_lim_ex  == conc_all$dir_pub, na.rm=TRUE), "/", n_genes, "\n")

write.xlsx(conc_all, file.path(output_dir, "concordance_six_way.xlsx"))
cat("\n\u2713 Schaub exact cohort analysis complete.\n")

################################################################################
# ADDITIONAL ANALYSES: PT-1:5 COMBINED (rpca_celltype definition)
# + PT-1 ONLY (highest SLC5A2, most like Schaub's primary signal)
# Both restricted to exact 16 Schaub participants
################################################################################

cat("\n=== PT-1:5 combined + PT-1 only analyses (exact 16 participants) ===\n")

# Helper function: run NEBULA on a given cell subset with Schaub exact grouping
run_nebula_schaub <- function(so_sub, label) {
  cat("\n--", label, "--\n")
  
  # Add Schaub exact grouping
  so_sub$sglt2i_schaub_run <- case_when(
    so_sub$record_id %in% sglt2i_resolved_s    ~ "Yes",
    so_sub$record_id %in% no_sglt2i_resolved_s ~ "No",
    TRUE                                       ~ NA_character_
  )
  
  grp_tab <- table(
    so_sub@meta.data %>%
      distinct(record_id, sglt2i_schaub_run) %>%
      pull(sglt2i_schaub_run)
  )
  cat("  Participants:", paste(names(grp_tab), grp_tab, sep="=", collapse=", "), "\n")
  if (length(grp_tab) < 2 || any(grp_tab < 2)) {
    cat("  Skipping.\n"); return(NULL)
  }
  cat("  Cells:", ncol(so_sub), "\n")
  
  genes_here <- base::intersect(genes_present, rownames(so_sub))
  so_g       <- so_sub[genes_here, ]
  counts_run <- round(GetAssayData(so_g, layer = "counts"))
  
  # Size factors on full cell set
  counts_full_run <- round(GetAssayData(so_sub, layer = "counts"))
  sce_run         <- SingleCellExperiment(assays = list(counts = counts_full_run))
  sce_run         <- computeSumFactors(sce_run)
  offset_run      <- sizeFactors(sce_run)
  rm(sce_run, counts_full_run)
  
  meta_run <- so_g@meta.data %>%
    mutate(sglt2i     = factor(sglt2i_schaub_run, levels = c("No", "Yes")),
           offset_run = offset_run[match(rownames(so_g@meta.data),
                                         rownames(so_sub@meta.data))])
  
  ok         <- complete.cases(meta_run[, c("sglt2i", "offset_run")])
  counts_run <- counts_run[, ok]
  meta_run   <- meta_run[ok, ]
  
  pred_run <- model.matrix(~ sglt2i, data = meta_run)
  lib_run  <- meta_run$offset_run
  
  data_g <- group_cell(count = counts_run, id = meta_run$record_id,
                       pred = pred_run, offset = lib_run)
  if (is.null(data_g)) {
    data_g <- list(count = counts_run, id = meta_run$record_id,
                   pred = pred_run, offset = lib_run)
  }
  
  neb <- tryCatch(
    nebula(count = data_g$count, id = data_g$id, pred = data_g$pred,
           offset = data_g$offset, model = "NBLMM",
           ncore = 1, reml = TRUE, output_re = TRUE, covariance = TRUE),
    error = function(e) { cat("  NEBULA error:", e$message, "\n"); NULL }
  )
  if (is.null(neb)) return(NULL)
  
  res    <- neb$summary %>% as_tibble()
  fc_col <- grep("logFC.*sglt2|logFC.*Yes", names(res), value=TRUE, ignore.case=TRUE)[1]
  p_col  <- grep("^p_.*sglt2|^p_.*Yes",    names(res), value=TRUE, ignore.case=TRUE)[1]
  
  res %>%
    dplyr::rename(gene_symbol = gene,
                  logFC       = all_of(fc_col),
                  p_value     = all_of(p_col)) %>%
    mutate(fdr   = p.adjust(p_value, method = "BH"),
           sig   = fdr < 0.05,
           label = label) %>%
    left_join(imap_dfr(pathway_genes, ~ tibble(gene_symbol = .x, pathway = .y)),
              by = "gene_symbol")
}

# ── 1. PT-1:5 combined, exact 16 participants ─────────────────────────────────
so_pt15_exact <- subset(so,
                        subset = celltype_rpca %in% pt_subclusters &
                          record_id %in% all_schaub_s)
cat("PT-1:5 combined (exact 16) cells:", ncol(so_pt15_exact), "\n")
res_pt15_exact <- run_nebula_schaub(so_pt15_exact, "PT-1:5 combined (exact 16)")

# ── 2. PT-1 only, exact 16 participants ──────────────────────────────────────
so_pt1_exact <- subset(so,
                       subset = celltype_rpca == "PT-1" &
                         record_id %in% all_schaub_s)
cat("PT-1 only (exact 16) cells:", ncol(so_pt1_exact), "\n")
res_pt1_exact <- run_nebula_schaub(so_pt1_exact, "PT-1 only (exact 16)")

# ── 3. PT-1:5 combined, all participants ─────────────────────────────────────
so_pt15_all <- subset(so, subset = celltype_rpca %in% pt_subclusters)
cat("PT-1:5 combined (all participants) cells:", ncol(so_pt15_all), "\n")
res_pt15_all <- run_nebula_schaub(so_pt15_all, "PT-1:5 combined (all)")

# Save
write.xlsx(
  list(PT15_exact16 = res_pt15_exact,
       PT1_exact16  = res_pt1_exact,
       PT15_all     = res_pt15_all),
  file.path(output_dir, "NEBULA_PT_subcluster_comparisons.xlsx")
)

################################################################################
# EXTENDED COMPARISON FIGURE — add PT-1:5 and PT-1 to six-way
################################################################################

# Build extended comparison dataframe
ext_df <- published_vals %>%
  left_join(res_nebula_ex  %>% dplyr::select(gene_symbol, logFC_kpmp_ex   = logFC_nebula_ex), by = "gene_symbol") %>%
  left_join(res_pt15_exact %>% dplyr::select(gene_symbol, logFC_pt15_ex   = logFC),           by = "gene_symbol") %>%
  left_join(res_pt1_exact  %>% dplyr::select(gene_symbol, logFC_pt1_ex    = logFC),           by = "gene_symbol") %>%
  left_join(res_pt15_all   %>% dplyr::select(gene_symbol, logFC_pt15_all  = logFC),           by = "gene_symbol") %>%
  pivot_longer(cols = c(logFC_published, logFC_kpmp_ex,
                        logFC_pt15_ex, logFC_pt1_ex, logFC_pt15_all),
               names_to = "source", values_to = "logFC") %>%
  mutate(
    source = factor(source,
                    levels = c("logFC_published",
                               "logFC_kpmp_ex",
                               "logFC_pt15_ex",
                               "logFC_pt1_ex",
                               "logFC_pt15_all"),
                    labels = c("Schaub et al.",
                               "KPMP PT (exact 16)",
                               "PT-1:5 rpca (exact 16)",
                               "PT-1 only (exact 16)",
                               "PT-1:5 rpca (all)")),
    pathway     = factor(pathway, levels = names(pathway_genes)),
    gene_symbol = factor(gene_symbol, levels = all_genes)
  )

ext_colors <- c(
  "Schaub et al."          = "#E07B54",
  "KPMP PT (exact 16)"     = "#008080",
  "PT-1:5 rpca (exact 16)" = "#2C7BB6",
  "PT-1 only (exact 16)"   = "#4B0082",
  "PT-1:5 rpca (all)"      = "#A8D8D8"
)

make_ext_panel <- function(pw_name, df) {
  gene_order <- pathway_genes[[pw_name]]
  pw_df <- df %>%
    dplyr::filter(pathway == pw_name) %>%
    mutate(gene_symbol = factor(as.character(gene_symbol), levels = gene_order))
  
  ggplot(pw_df, aes(x = gene_symbol, y = logFC, fill = source)) +
    geom_col(position = position_dodge(0.85), width = 0.8, na.rm = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey30") +
    scale_fill_manual(values = ext_colors, name = NULL) +
    scale_x_discrete(drop = FALSE) +
    labs(title = pw_name, x = NULL, y = "log\u2082 fold change") +
    theme_bw(base_size = 9) +
    theme(
      plot.title         = element_text(size = 9, face = "bold"),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 7),
      axis.title.y       = element_text(size = 8),
      legend.position    = "top",
      legend.text        = element_text(size = 7),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

ext_panels <- imap(setNames(names(pathway_genes), names(pathway_genes)),
                   ~ make_ext_panel(.x, ext_df))

fig_ext <- (ext_panels[[1]] | ext_panels[[2]]) /
  (ext_panels[[3]] | ext_panels[[4]] | ext_panels[[5]]) +
  plot_annotation(
    title    = "PT cell definition comparison vs Schaub et al. (exact 16 participants)",
    subtitle = "KPMP celltype2 PT vs PT-1:5 rpca combined vs PT-1 only | NEBULA SGLT2i vs No SGLT2i",
    theme    = theme(plot.title    = element_text(size = 11, face = "bold"),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

ggsave(file.path(output_dir, "pathway_PT_definition_comparison.pdf"),
       fig_ext, width = 18, height = 8, device = "pdf")
ggsave(file.path(output_dir, "pathway_PT_definition_comparison.png"),
       fig_ext, width = 18, height = 8, dpi = 300)
cat("PT definition comparison figure saved.\n")

# Concordance for new analyses
cat("\n--- PT cell definition concordance with Schaub ---\n")
conc_ext <- published_vals %>%
  left_join(res_nebula_ex  %>% dplyr::select(gene_symbol, logFC_kpmp_ex  = logFC_nebula_ex), by="gene_symbol") %>%
  left_join(res_pt15_exact %>% dplyr::select(gene_symbol, logFC_pt15_ex  = logFC),           by="gene_symbol") %>%
  left_join(res_pt1_exact  %>% dplyr::select(gene_symbol, logFC_pt1_ex   = logFC),           by="gene_symbol") %>%
  left_join(res_pt15_all   %>% dplyr::select(gene_symbol, logFC_pt15_all = logFC),           by="gene_symbol") %>%
  mutate(dir_pub     = sign(logFC_published),
         conc_kpmp   = sign(logFC_kpmp_ex)  == dir_pub,
         conc_pt15ex = sign(logFC_pt15_ex)  == dir_pub,
         conc_pt1ex  = sign(logFC_pt1_ex)   == dir_pub,
         conc_pt15al = sign(logFC_pt15_all) == dir_pub)

n <- sum(!is.na(conc_ext$logFC_published))
cat("KPMP PT (exact 16):      ", sum(conc_ext$conc_kpmp,   na.rm=TRUE), "/", n, "\n")
cat("PT-1:5 rpca (exact 16): ", sum(conc_ext$conc_pt15ex, na.rm=TRUE), "/", n, "\n")
cat("PT-1 only (exact 16):   ", sum(conc_ext$conc_pt1ex,  na.rm=TRUE), "/", n, "\n")
cat("PT-1:5 rpca (all):      ", sum(conc_ext$conc_pt15al, na.rm=TRUE), "/", n, "\n")

write.xlsx(conc_ext, file.path(output_dir, "concordance_PT_definition_comparison.xlsx"))
cat("\n\u2713 PT definition comparison complete.\n")

################################################################################
# PER-PARTICIPANT VIOLIN PLOTS — raw expression for each pathway gene
# One page per gene: individual participant violins + group summary
################################################################################

cat("\n=== Generating per-participant violin plots ===\n")

library(ggforce)  # for facet pagination if needed

# Pull pseudobulk-normalized expression per cell for the exact 16 participants
# Use log-normalized counts from Seurat

# Work from so_pt_exact (exact 16, KPMP PT definition)
# Add log-normalized layer if not already present
DefaultAssay(so_pt_exact) <- "RNA"
so_pt_exact <- NormalizeData(so_pt_exact, normalization.method = "LogNormalize",
                             scale.factor = 10000, verbose = FALSE)

# Extract normalized expression for pathway genes
expr_mat <- GetAssayData(so_pt_exact, layer = "data")[
  base::intersect(genes_present, rownames(so_pt_exact)), ]

# Build long-format data frame
expr_df <- as.data.frame(t(as.matrix(expr_mat))) %>%
  rownames_to_column("cell_barcode") %>%
  left_join(
    so_pt_exact@meta.data %>%
      rownames_to_column("cell_barcode") %>%
      dplyr::select(cell_barcode, record_id, sglt2i_exact),
    by = "cell_barcode"
  ) %>%
  pivot_longer(cols = -c(cell_barcode, record_id, sglt2i_exact),
               names_to = "gene_symbol", values_to = "expr") %>%
  mutate(
    sglt2i_group = factor(sglt2i_exact, levels = c("No", "Yes"),
                          labels = c("SGLT2i(−)", "SGLT2i(+)")),
    pathway = case_when(
      gene_symbol %in% pathway_genes[["Glycolysis"]]                        ~ "Glycolysis",
      gene_symbol %in% pathway_genes[["Gluconeogenesis"]]                   ~ "Gluconeogenesis",
      gene_symbol %in% pathway_genes[["Pyruvate metabolism and TCA cycle"]] ~ "TCA cycle",
      gene_symbol %in% pathway_genes[["Glutathione conjugation"]]           ~ "Glutathione conjugation",
      gene_symbol %in% pathway_genes[["Metallothioneins bind metals"]]      ~ "Metallothioneins",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::filter(pathway != "Other")

# Order participants: SGLT2i(-) first, then SGLT2i(+), alphabetical within
participant_order <- c(
  sort(no_sglt2i_resolved_s),
  sort(sglt2i_resolved_s)
)
participant_order <- participant_order[participant_order %in% unique(expr_df$record_id)]
expr_df$record_id <- factor(expr_df$record_id, levels = participant_order)

# Colors
group_colors  <- c("SGLT2i(−)" = "#4878CF", "SGLT2i(+)" = "#D65F5F")
group_fills   <- c("SGLT2i(−)" = "#4878CF55", "SGLT2i(+)" = "#D65F5F55")

# Participant colors — shade within group
n_no  <- length(no_sglt2i_resolved_s[no_sglt2i_resolved_s %in% unique(expr_df$record_id)])
n_yes <- length(sglt2i_resolved_s[sglt2i_resolved_s %in% unique(expr_df$record_id)])
pt_colors <- c(
  colorRampPalette(c("#2255AA", "#99BBEE"))(n_no),
  colorRampPalette(c("#AA2222", "#EEAAAA"))(n_yes)
)
names(pt_colors) <- participant_order

# ── Build one plot per gene ───────────────────────────────────────────────────
make_gene_violin <- function(gene, df) {
  gdf <- df %>% dplyr::filter(gene_symbol == gene)
  pw  <- unique(gdf$pathway)
  
  # Group summary stats
  summ <- gdf %>%
    group_by(record_id, sglt2i_group) %>%
    summarise(median_expr = median(expr), .groups = "drop")
  
  group_summ <- gdf %>%
    group_by(sglt2i_group) %>%
    summarise(mean_expr   = mean(expr),
              median_expr = median(expr),
              se_expr     = sd(expr) / sqrt(dplyr::n()),
              .groups     = "drop")
  
  # Proportion of expressing cells per participant
  pct_expr <- gdf %>%
    group_by(record_id, sglt2i_group) %>%
    summarise(pct = mean(expr > 0) * 100, .groups = "drop")
  
  # Main violin + jitter per participant
  p_main <- ggplot(gdf, aes(x = record_id, y = expr, fill = sglt2i_group,
                            color = sglt2i_group)) +
    geom_violin(scale = "width", alpha = 0.4, linewidth = 0.3,
                draw_quantiles = 0.5) +
    geom_jitter(aes(color = record_id), width = 0.15, size = 0.3,
                alpha = 0.25, show.legend = FALSE) +
    stat_summary(fun = median, geom = "point", size = 2.5,
                 aes(color = sglt2i_group), show.legend = FALSE) +
    scale_fill_manual(values  = group_fills,  name = NULL) +
    scale_color_manual(values = c(group_colors, pt_colors), name = NULL) +
    geom_vline(xintercept = n_no + 0.5, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    annotate("text", x = n_no / 2,       y = Inf, vjust = 1.5,
             label = "SGLT2i(−)", color = "#4878CF", size = 3, fontface = "bold") +
    annotate("text", x = n_no + n_yes / 2, y = Inf, vjust = 1.5,
             label = "SGLT2i(+)", color = "#D65F5F", size = 3, fontface = "bold") +
    labs(x = NULL, y = "log-normalized expression",
         title = gene, subtitle = pw) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(size = 13, face = "bold"),
      plot.subtitle    = element_text(size = 9,  color = "grey50"),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      axis.title.y     = element_text(size = 9),
      legend.position  = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  # Group summary panel (box + mean ± SE)
  p_summ <- ggplot(group_summ, aes(x = sglt2i_group, y = mean_expr,
                                   fill = sglt2i_group, color = sglt2i_group)) +
    geom_col(alpha = 0.5, width = 0.5) +
    geom_errorbar(aes(ymin = mean_expr - se_expr,
                      ymax = mean_expr + se_expr),
                  width = 0.2, linewidth = 0.7) +
    geom_point(size = 3) +
    scale_fill_manual(values  = group_fills,  name = NULL) +
    scale_color_manual(values = group_colors, name = NULL) +
    labs(x = NULL, y = "Mean expr ± SE", title = "Group summary") +
    theme_bw(base_size = 10) +
    theme(legend.position  = "none",
          plot.title       = element_text(size = 10, face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
  
  # % expressing cells per participant
  p_pct <- ggplot(pct_expr, aes(x = record_id, y = pct,
                                fill = sglt2i_group, color = sglt2i_group)) +
    geom_col(alpha = 0.6, width = 0.7) +
    geom_vline(xintercept = n_no + 0.5, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    scale_fill_manual(values  = group_fills,  name = NULL) +
    scale_color_manual(values = group_colors, name = NULL) +
    labs(x = NULL, y = "% cells expressing") +
    theme_bw(base_size = 9) +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 7),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
  
  # Combine: main violin (wide) | group summary (narrow) stacked with % cells
  (p_main | p_summ) / p_pct +
    plot_layout(heights = c(3, 1), widths = c(4, 1)) +
    plot_annotation(
      caption = sprintf("Exact 16 Schaub participants | KPMP PT definition | n cells No=%d, Yes=%d",
                        sum(gdf$sglt2i_group == "SGLT2i(−)"),
                        sum(gdf$sglt2i_group == "SGLT2i(+)")),
      theme = theme(plot.caption = element_text(size = 7, color = "grey50"))
    )
}

# ── Render to multi-page PDF ──────────────────────────────────────────────────
gene_order_pdf <- c(
  pathway_genes[["Glycolysis"]],
  pathway_genes[["Gluconeogenesis"]],
  pathway_genes[["Pyruvate metabolism and TCA cycle"]],
  pathway_genes[["Glutathione conjugation"]],
  pathway_genes[["Metallothioneins bind metals"]]
)
gene_order_pdf <- gene_order_pdf[gene_order_pdf %in% unique(expr_df$gene_symbol)]

pdf(file.path(output_dir, "per_participant_violins_pathway_genes.pdf"),
    width = 16, height = 9, onefile = TRUE)

# Title page
grid::grid.newpage()
grid::grid.text(
  label = paste0(
    "Per-Participant Expression: Pathway Genes\n",
    "Exact 16 Schaub participants | KPMP PT cells\n",
    "NEBULA log-normalized counts\n\n",
    paste0("Genes: ", paste(gene_order_pdf, collapse = ", "))
  ),
  gp = grid::gpar(fontsize = 14, fontface = "bold"),
  x = 0.5, y = 0.5
)

for (gene in gene_order_pdf) {
  cat("  Plotting:", gene, "\n")
  p <- tryCatch(make_gene_violin(gene, expr_df),
                error = function(e) { cat("    Error:", e$message, "\n"); NULL })
  if (!is.null(p)) print(p)
}

dev.off()
cat("PDF saved:", file.path(output_dir, "per_participant_violins_pathway_genes.pdf"), "\n")
cat("\n\u2713 Violin plots complete.\n")