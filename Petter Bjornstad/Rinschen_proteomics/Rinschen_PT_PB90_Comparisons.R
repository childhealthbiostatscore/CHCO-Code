################################################################################
# PT PATHWAY SGLT2i ANALYSIS — PB_90_RPCAFix Seurat object
# Exact Schaub et al. 2023 (JCI) 16 T2D participants
# NEBULA + limma pseudobulk + per-participant violins
################################################################################

library(Seurat)
library(tidyverse)
library(tibble)
library(patchwork)
library(nebula)
library(SingleCellExperiment)
library(scran)
library(limma)
library(edgeR)
library(openxlsx)
library(conflicted)
library(gtsummary)

conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::rename)
conflicted::conflicts_prefer(dplyr::count)
conflicted::conflicts_prefer(base::intersect)
conflicted::conflicts_prefer(base::setdiff)
conflicted::conflicts_prefer(Seurat::Assays)

# ── Paths ─────────────────────────────────────────────────────────────────────
so_path    <- "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/scRNA/data_raw/PB_90_RPCAFix_Metadata_LR_RM_kpmpV1labelled.rds"
schaub_csv <- "C:/Users/netio/Downloads/jci2023_schaub_ids.csv"
output_dir <- "C:/Users/netio/Documents/UofW/proteomics_vs_PET/PB90_analysis/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# LOAD AND EXPLORE NEW SEURAT OBJECT
################################################################################

cat("Loading Seurat object...\n")
so <- readRDS(so_path)
cat("Object loaded:", ncol(so), "cells,", nrow(so), "genes\n")
cat("\nDimensions:", dim(so), "\n")

cat("\n--- Meta.data columns ---\n")
print(names(so@meta.data))

cat("\n--- First few rows of metadata ---\n")
print(head(so@meta.data[, 1:min(10, ncol(so@meta.data))]))

# Find cell type and participant ID columns
cat("\n--- Unique values for likely cell type columns ---\n")
possible_ct_cols <- names(so@meta.data)[
  str_detect(names(so@meta.data),
             regex("celltype|cell_type|cluster|ident|annotation|kpmp|rpca",
                   ignore_case = TRUE))
]
cat("Candidate celltype columns:", paste(possible_ct_cols, collapse = ", "), "\n")
for (col in possible_ct_cols) {
  cat("\n", col, ":\n")
  print(table(so@meta.data[[col]]))
}

cat("\n--- Unique values for likely participant ID columns ---\n")
possible_id_cols <- names(so@meta.data)[
  str_detect(names(so@meta.data),
             regex("record|participant|subject|id|sample|donor",
                   ignore_case = TRUE))
]
cat("Candidate ID columns:", paste(possible_id_cols, collapse = ", "), "\n")
for (col in possible_id_cols) {
  n_unique <- length(unique(so@meta.data[[col]]))
  if (n_unique < 100) {
    cat("\n", col, "(", n_unique, "unique):\n")
    print(sort(unique(so@meta.data[[col]])))
  }
}

cat("\n--- SGLT2i-related columns ---\n")
sglt2_cols <- names(so@meta.data)[
  str_detect(names(so@meta.data), regex("sglt2|sglti2", ignore_case = TRUE))
]
cat("SGLT2i columns:", paste(sglt2_cols, collapse = ", "), "\n")
for (col in sglt2_cols) {
  cat("\n", col, ":\n")
  print(table(so@meta.data[[col]], useNA = "ifany"))
}

cat("\n--- Available assays ---\n")
print(Assays(so))
cat("Default assay:", DefaultAssay(so), "\n")

################################################################################
# STOP HERE — run the above first to identify the correct column names,
# then update the variables below and continue
################################################################################

# ── UPDATE THESE based on the output above ────────────────────────────────────
CT_COL     <- "celltype_rpca"   # PT-1 through PT-5
CT_KPMP    <- "KPMP_celltype"   # PT-S1/S2, aPT, PT-S3
ID_COL     <- "record_id"       # RH-xx-T / IT_xx format
SGLT2_COL  <- "sglt2i_ever"     # Yes/No (treat "." as NA)

cat("\n=== Using columns: CT =", CT_COL, "| ID =", ID_COL,
    "| SGLT2i =", SGLT2_COL, "===\n")

# Fix "." coded as missing in sglt2i_ever
so@meta.data[[SGLT2_COL]] <- ifelse(
  so@meta.data[[SGLT2_COL]] == ".", NA_character_,
  so@meta.data[[SGLT2_COL]]
)
cat("sglt2i_ever after NA fix:\n")
print(table(so@meta.data[[SGLT2_COL]], useNA = "ifany"))

# Confirm columns exist
stopifnot(
  CT_COL   %in% names(so@meta.data),
  ID_COL   %in% names(so@meta.data),
  SGLT2_COL %in% names(so@meta.data)
)

################################################################################
# SCHAUB EXACT 16 PARTICIPANTS
################################################################################

schaub_meta <- read.csv(schaub_csv, stringsAsFactors = FALSE) %>%
  distinct(ID, T2D_HC, SGLT2i, SGLT2i_duration_in_exposed) %>%
  mutate(group_schaub = case_when(
    T2D_HC == "HC"              ~ "HC",
    SGLT2i == "YES"             ~ "SGLT2i",
    T2D_HC %in% c("T2D","T2Di") & SGLT2i == "NO" ~ "No SGLT2i",
    TRUE ~ NA_character_
  ))

schaub_sglt2i_ids    <- schaub_meta %>%
  dplyr::filter(SGLT2i == "YES") %>% pull(ID)
schaub_no_sglt2i_ids <- schaub_meta %>%
  dplyr::filter(T2D_HC %in% c("T2D","T2Di"), SGLT2i == "NO") %>% pull(ID)

cat("Schaub SGLT2i(+):", paste(sort(schaub_sglt2i_ids), collapse=", "), "\n")
cat("Schaub SGLT2i(-):", paste(sort(schaub_no_sglt2i_ids), collapse=", "), "\n")

# Resolve co-enrolled aliases
coenroll_alias <- c("IT_07"="RH-59-T","IT_08"="RH-60-T",
                    "IT_09"="RH-65-T","IT_10"="RH-66-T")

resolve_id <- function(ids, so_obj) {
  all_ids <- unique(so_obj@meta.data[[ID_COL]])
  sapply(ids, function(id) {
    if (id %in% all_ids) return(id)
    alias <- coenroll_alias[id]
    if (!is.na(alias) && alias %in% all_ids) return(alias)
    return(NA_character_)
  })
}

sglt2i_resolved    <- na.omit(resolve_id(schaub_sglt2i_ids,    so))
no_sglt2i_resolved <- na.omit(resolve_id(schaub_no_sglt2i_ids, so))
all_schaub_ids     <- unique(c(sglt2i_resolved, no_sglt2i_resolved))

cat("\nResolved in new object:\n")
cat("SGLT2i(+):", paste(sort(sglt2i_resolved),    collapse=", "), "\n")
cat("SGLT2i(-):", paste(sort(no_sglt2i_resolved), collapse=", "), "\n")
cat("Total:", length(all_schaub_ids), "/ 16\n")

# Any missing?
missing_ids <- base::setdiff(
  c(schaub_sglt2i_ids, schaub_no_sglt2i_ids),
  c(sglt2i_resolved, no_sglt2i_resolved)
)
if (length(missing_ids) > 0)
  cat("WARNING — not found:", paste(missing_ids, collapse=", "), "\n")

################################################################################
# PATHWAY GENE SETS (from Schaub Fig 4C)
################################################################################

pathway_genes <- list(
  "Glycolysis" = c(
    "PKLR","PFKFB3","PFKL","ALDOC","HK2","ENO2",
    "PGK1","PGAM1","TPI1","GAPDH"
  ),
  "Gluconeogenesis" = c(
    "SLC25A10","GOT2","GOT1","FBP1","SLC25A11","PCK1","MDH1"
  ),
  "Pyruvate metabolism and TCA cycle" = c(
    "SDHB","SUCLG1","PDK2","ACO2","IDH3G",
    "SUCLA2","HAGH","PDHB","LDHA"
  ),
  "Glutathione conjugation" = c(
    "CNDP2","GSTM4","GSTT2B","GSTO1","GGCT","GSTM3","AKR1A1"
  ),
  "Metallothioneins bind metals" = c(
    "MT1G","MT1X","MT1H","MT2A"
  )
)
all_genes <- unlist(pathway_genes, use.names = FALSE)

# Published values from Schaub Fig 4C (corrected from paper)
published_vals <- tribble(
  ~gene_symbol,  ~pathway,                                ~logFC_published,
  "PKLR",        "Glycolysis",                            -0.02,
  "PFKFB3",      "Glycolysis",                            -0.02,
  "PFKL",        "Glycolysis",                            -0.02,
  "ALDOC",       "Glycolysis",                            -0.02,
  "HK2",         "Glycolysis",                            -0.03,
  "ENO2",        "Glycolysis",                            -0.02,
  "PGK1",        "Glycolysis",                             0.05,
  "PGAM1",       "Glycolysis",                             0.05,
  "TPI1",        "Glycolysis",                             0.13,
  "GAPDH",       "Glycolysis",                            -0.04,
  "SLC25A10",    "Gluconeogenesis",                        0.00,
  "GOT2",        "Gluconeogenesis",                        0.05,
  "GOT1",        "Gluconeogenesis",                        0.05,
  "FBP1",        "Gluconeogenesis",                       -0.05,
  "SLC25A11",    "Gluconeogenesis",                       -0.10,
  "PCK1",        "Gluconeogenesis",                       -0.15,
  "MDH1",        "Gluconeogenesis",                       -0.15,
  "SDHB",        "Pyruvate metabolism and TCA cycle",     -0.05,
  "SUCLG1",      "Pyruvate metabolism and TCA cycle",     -0.05,
  "PDK2",        "Pyruvate metabolism and TCA cycle",     -0.05,
  "ACO2",        "Pyruvate metabolism and TCA cycle",     -0.05,
  "IDH3G",       "Pyruvate metabolism and TCA cycle",     -0.05,
  "SUCLA2",      "Pyruvate metabolism and TCA cycle",     -0.08,
  "HAGH",        "Pyruvate metabolism and TCA cycle",     -0.05,
  "PDHB",        "Pyruvate metabolism and TCA cycle",     -0.08,
  "LDHA",        "Pyruvate metabolism and TCA cycle",     -0.10,
  "CNDP2",       "Glutathione conjugation",               -0.02,
  "GSTM4",       "Glutathione conjugation",               -0.03,
  "GSTT2B",      "Glutathione conjugation",               -0.03,
  "GSTO1",       "Glutathione conjugation",               -0.03,
  "GGCT",        "Glutathione conjugation",                0.00,
  "GSTM3",       "Glutathione conjugation",               -0.05,
  "AKR1A1",      "Glutathione conjugation",                0.08,
  "MT1G",        "Metallothioneins bind metals",           0.08,
  "MT1X",        "Metallothioneins bind metals",           0.08,
  "MT1H",        "Metallothioneins bind metals",           0.08,
  "MT2A",        "Metallothioneins bind metals",           0.08
)

################################################################################
# PT CELL SUBSETS — using new object's cell type column
################################################################################

# Identify PT cells — check what labels exist
cat("\nCell type labels in new object:\n")
print(sort(unique(so@meta.data[[CT_COL]])))

# PT-1:5 rpca subclusters (celltype_rpca)
PT_RPCA_LABELS <- c("PT-1","PT-2","PT-3","PT-4","PT-5")

# KPMP celltype labels (KPMP_celltype column in this object)
PT_KPMP_LABELS <- c("PT-S1/S2","aPT","PT-S3")

# Filter to what actually exists
pt_rpca_present <- base::intersect(PT_RPCA_LABELS,
                                   unique(so@meta.data[[CT_COL]]))
pt_kpmp_present <- base::intersect(PT_KPMP_LABELS,
                                   unique(so@meta.data[[CT_KPMP]]))

cat("PT rpca labels found:", paste(pt_rpca_present, collapse=", "), "\n")
cat("PT KPMP labels found:", paste(pt_kpmp_present, collapse=", "), "\n")

################################################################################
# HELPER: run NEBULA on a Seurat subset
################################################################################

run_nebula_subset <- function(so_sub, label,
                              sglt2i_yes_ids, sglt2i_no_ids,
                              genes_vec, id_col, sglt2i_col) {
  cat("\n--", label, "--\n")
  cat("  Cells:", ncol(so_sub), "\n")
  
  # Add group assignment
  so_sub$sglt2i_run <- case_when(
    so_sub@meta.data[[id_col]] %in% sglt2i_yes_ids ~ "Yes",
    so_sub@meta.data[[id_col]] %in% sglt2i_no_ids  ~ "No",
    TRUE ~ NA_character_
  )
  
  grp_tab <- table(
    so_sub@meta.data %>%
      distinct(across(all_of(id_col)), sglt2i_run) %>%
      pull(sglt2i_run)
  )
  cat("  Participants:", paste(names(grp_tab), grp_tab, sep="=", collapse=", "), "\n")
  
  if (length(grp_tab) < 2 || any(grp_tab < 2)) {
    cat("  Skipping — insufficient groups.\n"); return(NULL)
  }
  
  genes_here <- base::intersect(genes_vec, rownames(so_sub))
  cat("  Genes found:", length(genes_here), "/", length(genes_vec), "\n")
  so_g <- so_sub[genes_here, ]
  counts_sub <- round(GetAssayData(so_g, layer = "counts"))
  
  # Size factors
  counts_full <- round(GetAssayData(so_sub, layer = "counts"))
  sce_tmp     <- SingleCellExperiment(assays = list(counts = counts_full))
  sce_tmp     <- computeSumFactors(sce_tmp)
  offset_vec  <- sizeFactors(sce_tmp)
  rm(sce_tmp, counts_full)
  
  meta_sub <- so_g@meta.data %>%
    mutate(sglt2i    = factor(sglt2i_run, levels = c("No","Yes")),
           offset_sf = offset_vec[match(rownames(so_g@meta.data),
                                        rownames(so_sub@meta.data))])
  
  ok         <- complete.cases(meta_sub[, c("sglt2i","offset_sf")])
  counts_sub <- counts_sub[, ok]
  meta_sub   <- meta_sub[ok, ]
  
  pred_sub <- model.matrix(~ sglt2i, data = meta_sub)
  lib_sub  <- meta_sub$offset_sf
  
  data_g <- group_cell(count = counts_sub,
                       id    = meta_sub[[id_col]],
                       pred  = pred_sub, offset = lib_sub)
  if (is.null(data_g)) {
    data_g <- list(count = counts_sub, id = meta_sub[[id_col]],
                   pred  = pred_sub,   offset = lib_sub)
  }
  
  neb <- tryCatch(
    nebula(count = data_g$count, id = data_g$id, pred = data_g$pred,
           offset = data_g$offset, model = "NBLMM",
           ncore = 1, reml = TRUE, output_re = TRUE, covariance = TRUE),
    error = function(e) { cat("  NEBULA error:", e$message, "\n"); NULL }
  )
  if (is.null(neb)) return(NULL)
  
  res    <- neb$summary %>% as_tibble()
  fc_col <- grep("logFC.*sglt2|logFC.*Yes", names(res), value=TRUE,
                 ignore.case=TRUE)[1]
  p_col  <- grep("^p_.*sglt2|^p_.*Yes", names(res), value=TRUE,
                 ignore.case=TRUE)[1]
  
  res %>%
    dplyr::rename(gene_symbol = gene,
                  logFC       = all_of(fc_col),
                  p_value     = all_of(p_col)) %>%
    mutate(fdr   = p.adjust(p_value, method = "BH"),
           sig   = fdr < 0.05,
           label = label) %>%
    left_join(imap_dfr(pathway_genes, ~ tibble(gene_symbol=.x, pathway=.y)),
              by = "gene_symbol")
}

################################################################################
# HELPER: run limma pseudobulk
################################################################################

run_limma_subset <- function(so_sub, label,
                             sglt2i_yes_ids, sglt2i_no_ids,
                             genes_vec, id_col) {
  cat("\n-- Limma:", label, "--\n")
  
  so_sub$sglt2i_run <- case_when(
    so_sub@meta.data[[id_col]] %in% sglt2i_yes_ids ~ "Yes",
    so_sub@meta.data[[id_col]] %in% sglt2i_no_ids  ~ "No",
    TRUE ~ NA_character_
  )
  
  genes_here   <- base::intersect(genes_vec, rownames(so_sub))
  counts_cells <- round(GetAssayData(so_sub[genes_here,], layer = "counts"))
  meta_cells   <- so_sub@meta.data %>%
    mutate(sglt2i = factor(sglt2i_run, levels = c("No","Yes"))) %>%
    dplyr::select(all_of(id_col), sglt2i)
  
  pids <- unique(meta_cells[[id_col]])
  pb   <- sapply(pids, function(pid) {
    cells <- which(meta_cells[[id_col]] == pid)
    Matrix::rowSums(counts_cells[, cells, drop=FALSE])
  })
  colnames(pb) <- pids
  
  pb_meta <- meta_cells %>%
    distinct(across(all_of(id_col)), .keep_all=TRUE) %>%
    dplyr::filter(!is.na(sglt2i)) %>%
    arrange(match(.data[[id_col]], pids))
  pb <- pb[, pb_meta[[id_col]]]
  
  cat("  Pseudobulk:", nrow(pb), "genes x", ncol(pb), "participants\n")
  cat("  SGLT2i:", table(pb_meta$sglt2i), "\n")
  
  dge    <- DGEList(counts=pb, group=pb_meta$sglt2i)
  dge    <- calcNormFactors(dge, method="TMM")
  design <- model.matrix(~ sglt2i, data=pb_meta)
  v      <- voom(dge, design, plot=FALSE)
  fit    <- eBayes(lmFit(v, design))
  
  topTable(fit, coef="sglt2iYes", number=Inf, sort.by="none") %>%
    rownames_to_column("gene_symbol") %>% as_tibble() %>%
    dplyr::rename(logFC_limma = logFC, fdr_limma = adj.P.Val) %>%
    mutate(sig_limma = fdr_limma < 0.05, label = label) %>%
    left_join(imap_dfr(pathway_genes, ~ tibble(gene_symbol=.x, pathway=.y)),
              by="gene_symbol")
}

################################################################################
# GENES PRESENT CHECK
################################################################################

genes_present <- base::intersect(all_genes, rownames(so))
cat("\nGenes present in new object:", length(genes_present), "/", length(all_genes), "\n")
missing_genes <- base::setdiff(all_genes, rownames(so))
if (length(missing_genes) > 0)
  cat("Missing genes:", paste(missing_genes, collapse=", "), "\n")

################################################################################
# RUN ALL ANALYSES
################################################################################

results_list <- list()

# ── 1. KPMP PT definition (PT-S1/S2, aPT, PT-S3), exact 16 ──────────────────
if (length(pt_kpmp_present) > 0) {
  so_pt_kpmp_ex <- subset(so,
                          cells = colnames(so)[
                            so@meta.data[[CT_KPMP]] %in% pt_kpmp_present &
                              so@meta.data[[ID_COL]]  %in% all_schaub_ids
                          ])
  results_list[["NEBULA_KPMP_exact16"]] <- run_nebula_subset(
    so_pt_kpmp_ex, "KPMP PT (exact 16)",
    sglt2i_resolved, no_sglt2i_resolved,
    genes_present, ID_COL, SGLT2_COL
  )
  results_list[["Limma_KPMP_exact16"]] <- run_limma_subset(
    so_pt_kpmp_ex, "KPMP PT limma (exact 16)",
    sglt2i_resolved, no_sglt2i_resolved,
    genes_present, ID_COL
  )
}

# ── 2. PT-1:5 rpca combined, exact 16 ────────────────────────────────────────
if (length(pt_rpca_present) > 0) {
  so_pt15_ex <- subset(so,
                       cells = colnames(so)[
                         so@meta.data[[CT_COL]] %in% pt_rpca_present &
                           so@meta.data[[ID_COL]] %in% all_schaub_ids
                       ])
  results_list[["NEBULA_PT15_exact16"]] <- run_nebula_subset(
    so_pt15_ex, "PT-1:5 rpca (exact 16)",
    sglt2i_resolved, no_sglt2i_resolved,
    genes_present, ID_COL, SGLT2_COL
  )
  results_list[["Limma_PT15_exact16"]] <- run_limma_subset(
    so_pt15_ex, "PT-1:5 rpca limma (exact 16)",
    sglt2i_resolved, no_sglt2i_resolved,
    genes_present, ID_COL
  )
  
  # ── 3. PT-1 only, exact 16 ─────────────────────────────────────────────────
  if ("PT-1" %in% pt_rpca_present) {
    so_pt1_ex <- subset(so,
                        cells = colnames(so)[
                          so@meta.data[[CT_COL]] == "PT-1" &
                            so@meta.data[[ID_COL]] %in% all_schaub_ids
                        ])
    results_list[["NEBULA_PT1_exact16"]] <- run_nebula_subset(
      so_pt1_ex, "PT-1 only (exact 16)",
      sglt2i_resolved, no_sglt2i_resolved,
      genes_present, ID_COL, SGLT2_COL
    )
  }
  
  # ── 4. PT-1:5 rpca combined, ALL participants ─────────────────────────────
  so_pt15_all <- subset(so,
                        cells = colnames(so)[so@meta.data[[CT_COL]] %in% pt_rpca_present])
  results_list[["NEBULA_PT15_all"]] <- run_nebula_subset(
    so_pt15_all, "PT-1:5 rpca (all)",
    unique(so@meta.data[[ID_COL]][so@meta.data[[SGLT2_COL]] == "Yes" &
                                    !is.na(so@meta.data[[SGLT2_COL]])]),
    unique(so@meta.data[[ID_COL]][so@meta.data[[SGLT2_COL]] == "No"  &
                                    !is.na(so@meta.data[[SGLT2_COL]])]),
    genes_present, ID_COL, SGLT2_COL
  )
}

# Save all DE results
write.xlsx(compact(results_list),
           file.path(output_dir, "PB90_NEBULA_limma_pathway_results.xlsx"))
cat("\nResults saved.\n")

################################################################################
# COMPARISON FIGURE vs Schaub
################################################################################

# Pull logFC for each analysis
get_lfc <- function(res, gene_col = "gene_symbol", lfc_col = "logFC") {
  if (is.null(res)) return(tibble(gene_symbol = character(), logFC = numeric()))
  res %>% dplyr::select(gene_symbol = all_of(gene_col), logFC = all_of(lfc_col))
}

comp_df <- published_vals %>%
  left_join(get_lfc(results_list[["NEBULA_KPMP_exact16"]]) %>%
              dplyr::rename(logFC_kpmp_ex = logFC),    by="gene_symbol") %>%
  left_join(get_lfc(results_list[["NEBULA_PT15_exact16"]]) %>%
              dplyr::rename(logFC_pt15_ex = logFC),    by="gene_symbol") %>%
  left_join(get_lfc(results_list[["NEBULA_PT1_exact16"]]) %>%
              dplyr::rename(logFC_pt1_ex  = logFC),    by="gene_symbol") %>%
  left_join(get_lfc(results_list[["NEBULA_PT15_all"]]) %>%
              dplyr::rename(logFC_pt15_all = logFC),   by="gene_symbol") %>%
  left_join(get_lfc(results_list[["Limma_KPMP_exact16"]], lfc_col="logFC_limma") %>%
              dplyr::rename(logFC_limma_ex = logFC),   by="gene_symbol") %>%
  pivot_longer(cols = starts_with("logFC_") | matches("logFC_published"),
               names_to = "source", values_to = "logFC") %>%
  mutate(
    source = factor(source,
                    levels = c("logFC_published","logFC_kpmp_ex","logFC_pt15_ex",
                               "logFC_pt1_ex","logFC_pt15_all","logFC_limma_ex"),
                    labels = c("Schaub et al.","KPMP PT (exact 16)",
                               "PT-1:5 rpca (exact 16)","PT-1 only (exact 16)",
                               "PT-1:5 rpca (all)","Limma KPMP (exact 16)")),
    pathway     = factor(pathway, levels = names(pathway_genes)),
    gene_symbol = factor(gene_symbol, levels = all_genes)
  )

comp_colors <- c(
  "Schaub et al."          = "#E07B54",
  "KPMP PT (exact 16)"     = "#008080",
  "PT-1:5 rpca (exact 16)" = "#2C7BB6",
  "PT-1 only (exact 16)"   = "#4B0082",
  "PT-1:5 rpca (all)"      = "#A8D8D8",
  "Limma KPMP (exact 16)"  = "#C4A8E0"
)

make_comp_panel <- function(pw_name, df) {
  pw_df <- df %>%
    dplyr::filter(pathway == pw_name) %>%
    mutate(gene_symbol = factor(as.character(gene_symbol),
                                levels = pathway_genes[[pw_name]]))
  ggplot(pw_df, aes(x = gene_symbol, y = logFC, fill = source)) +
    geom_col(position = position_dodge(0.85), width = 0.8, na.rm = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey30") +
    scale_fill_manual(values = comp_colors, name = NULL) +
    scale_x_discrete(drop = FALSE) +
    labs(title = pw_name, x = NULL, y = "log\u2082 fold change") +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.title.y = element_text(size = 8),
          legend.position = "top", legend.text = element_text(size = 7),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
}

panels <- imap(setNames(names(pathway_genes), names(pathway_genes)),
               ~ make_comp_panel(.x, comp_df))

fig <- (panels[[1]] | panels[[2]]) /
  (panels[[3]] | panels[[4]] | panels[[5]]) +
  plot_annotation(
    title    = "PB90 object: Pathway comparison vs Schaub et al.",
    subtitle = "NEBULA SGLT2i vs No SGLT2i | exact 16 Schaub participants",
    theme    = theme(plot.title    = element_text(size = 11, face = "bold"),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

ggsave(file.path(output_dir, "PB90_pathway_comparison.pdf"),
       fig, width = 18, height = 8, device = "pdf")
ggsave(file.path(output_dir, "PB90_pathway_comparison.png"),
       fig, width = 18, height = 8, dpi = 300)
cat("Comparison figure saved.\n")

################################################################################
# PER-PARTICIPANT VIOLIN PLOTS
################################################################################

cat("\n=== Generating per-participant violin plots ===\n")

# Use KPMP PT exact 16 if available, else PT-1:5
so_vln <- if (length(pt_kpmp_present) > 0) so_pt_kpmp_ex else so_pt15_ex
so_vln$sglt2i_group_vln <- case_when(
  so_vln@meta.data[[ID_COL]] %in% sglt2i_resolved    ~ "Yes",
  so_vln@meta.data[[ID_COL]] %in% no_sglt2i_resolved ~ "No",
  TRUE ~ NA_character_
)
so_vln <- subset(so_vln, !is.na(sglt2i_group_vln))

DefaultAssay(so_vln) <- "RNA"
so_vln <- NormalizeData(so_vln, normalization.method="LogNormalize",
                        scale.factor=10000, verbose=FALSE)

expr_mat <- GetAssayData(so_vln, layer="data")[
  base::intersect(genes_present, rownames(so_vln)), ]

expr_df <- as.data.frame(t(as.matrix(expr_mat))) %>%
  rownames_to_column("cell_barcode") %>%
  left_join(
    so_vln@meta.data %>%
      rownames_to_column("cell_barcode") %>%
      dplyr::select(cell_barcode,
                    record_id = all_of(ID_COL),
                    sglt2i_grp = sglt2i_group_vln),
    by = "cell_barcode"
  ) %>%
  pivot_longer(cols = -c(cell_barcode, record_id, sglt2i_grp),
               names_to = "gene_symbol", values_to = "expr") %>%
  mutate(
    sglt2i_group = factor(sglt2i_grp, levels=c("No","Yes"),
                          labels=c("SGLT2i(\u2212)","SGLT2i(+)")),
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

# Participant order
pt_order <- c(sort(no_sglt2i_resolved), sort(sglt2i_resolved))
pt_order  <- pt_order[pt_order %in% unique(expr_df$record_id)]
expr_df$record_id <- factor(expr_df$record_id, levels = pt_order)

n_no  <- sum(pt_order %in% no_sglt2i_resolved)
n_yes <- sum(pt_order %in% sglt2i_resolved)

group_colors <- c("SGLT2i(\u2212)"="#4878CF","SGLT2i(+)"="#D65F5F")
group_fills  <- c("SGLT2i(\u2212)"="#4878CF55","SGLT2i(+)"="#D65F5F55")
pt_colors    <- c(colorRampPalette(c("#2255AA","#99BBEE"))(n_no),
                  colorRampPalette(c("#AA2222","#EEAAAA"))(n_yes))
names(pt_colors) <- pt_order

make_gene_violin <- function(gene, df) {
  gdf  <- df %>% dplyr::filter(gene_symbol == gene)
  pw   <- unique(gdf$pathway)
  
  group_summ <- gdf %>%
    group_by(sglt2i_group) %>%
    summarise(mean_expr = mean(expr),
              se_expr   = sd(expr)/sqrt(dplyr::n()),
              .groups   = "drop")
  
  pct_expr <- gdf %>%
    group_by(record_id, sglt2i_group) %>%
    summarise(pct = mean(expr > 0)*100, .groups="drop")
  
  p_main <- ggplot(gdf, aes(x=record_id, y=expr,
                            fill=sglt2i_group, color=sglt2i_group)) +
    geom_violin(scale="width", alpha=0.4, linewidth=0.3, draw_quantiles=0.5) +
    geom_jitter(aes(color=record_id), width=0.15, size=0.3,
                alpha=0.25, show.legend=FALSE) +
    stat_summary(fun=median, geom="point", size=2.5,
                 aes(color=sglt2i_group), show.legend=FALSE) +
    scale_fill_manual(values=group_fills,  name=NULL) +
    scale_color_manual(values=c(group_colors, pt_colors), name=NULL) +
    geom_vline(xintercept=n_no+0.5, linetype="dashed",
               color="grey50", linewidth=0.5) +
    annotate("text", x=n_no/2,         y=Inf, vjust=1.5,
             label="SGLT2i(\u2212)", color="#4878CF", size=3, fontface="bold") +
    annotate("text", x=n_no+n_yes/2,   y=Inf, vjust=1.5,
             label="SGLT2i(+)",     color="#D65F5F", size=3, fontface="bold") +
    labs(x=NULL, y="log-normalized expression", title=gene, subtitle=pw) +
    theme_bw(base_size=10) +
    theme(plot.title=element_text(size=13, face="bold"),
          plot.subtitle=element_text(size=9, color="grey50"),
          axis.text.x=element_text(angle=45, hjust=1, size=8),
          legend.position="none",
          panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank())
  
  p_summ <- ggplot(group_summ,
                   aes(x=sglt2i_group, y=mean_expr,
                       fill=sglt2i_group, color=sglt2i_group)) +
    geom_col(alpha=0.5, width=0.5) +
    geom_errorbar(aes(ymin=mean_expr-se_expr, ymax=mean_expr+se_expr),
                  width=0.2, linewidth=0.7) +
    geom_point(size=3) +
    scale_fill_manual(values=group_fills,  name=NULL) +
    scale_color_manual(values=group_colors, name=NULL) +
    labs(x=NULL, y="Mean ± SE", title="Group summary") +
    theme_bw(base_size=10) +
    theme(legend.position="none", plot.title=element_text(size=10,face="bold"),
          panel.grid.minor=element_blank(), panel.grid.major.x=element_blank())
  
  p_pct <- ggplot(pct_expr, aes(x=record_id, y=pct,
                                fill=sglt2i_group, color=sglt2i_group)) +
    geom_col(alpha=0.6, width=0.7) +
    geom_vline(xintercept=n_no+0.5, linetype="dashed",
               color="grey50", linewidth=0.5) +
    scale_fill_manual(values=group_fills,  name=NULL) +
    scale_color_manual(values=group_colors, name=NULL) +
    labs(x=NULL, y="% cells expressing") +
    theme_bw(base_size=9) +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=7),
          legend.position="none",
          panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank())
  
  (p_main | p_summ) / p_pct +
    plot_layout(heights=c(3,1), widths=c(4,1)) +
    plot_annotation(
      caption = sprintf("PB90 object | exact 16 Schaub participants | n cells No=%d, Yes=%d",
                        sum(gdf$sglt2i_group=="SGLT2i(\u2212)"),
                        sum(gdf$sglt2i_group=="SGLT2i(+)")),
      theme = theme(plot.caption=element_text(size=7, color="grey50"))
    )
}

gene_order_pdf <- unlist(pathway_genes, use.names=FALSE)
gene_order_pdf <- gene_order_pdf[gene_order_pdf %in% unique(expr_df$gene_symbol)]

pdf(file.path(output_dir, "PB90_per_participant_violins.pdf"),
    width=16, height=9, onefile=TRUE)
grid::grid.newpage()
grid::grid.text(
  paste0("PB90 Object: Per-Participant Expression\n",
         "Exact 16 Schaub participants\n",
         paste(gene_order_pdf, collapse=", ")),
  gp=grid::gpar(fontsize=14, fontface="bold"), x=0.5, y=0.5
)
for (gene in gene_order_pdf) {
  cat("  Plotting:", gene, "\n")
  p <- tryCatch(make_gene_violin(gene, expr_df),
                error=function(e){cat("Error:",e$message,"\n");NULL})
  if (!is.null(p)) print(p)
}
dev.off()
cat("Violin PDF saved.\n")

# Concordance summary
cat("\n--- Concordance with Schaub et al. ---\n")
conc <- published_vals %>%
  left_join(results_list[["NEBULA_KPMP_exact16"]]  %>% dplyr::select(gene_symbol, logFC_kpmp  = logFC), by="gene_symbol") %>%
  left_join(results_list[["NEBULA_PT15_exact16"]]  %>% dplyr::select(gene_symbol, logFC_pt15  = logFC), by="gene_symbol") %>%
  left_join(results_list[["NEBULA_PT1_exact16"]]   %>% dplyr::select(gene_symbol, logFC_pt1   = logFC), by="gene_symbol") %>%
  mutate(dir_pub  = sign(logFC_published),
         conc_kpmp = sign(logFC_kpmp)  == dir_pub,
         conc_pt15 = sign(logFC_pt15)  == dir_pub,
         conc_pt1  = sign(logFC_pt1)   == dir_pub)
n <- sum(!is.na(conc$logFC_published))
cat("KPMP PT (exact 16):      ", sum(conc$conc_kpmp, na.rm=TRUE), "/", n, "\n")
cat("PT-1:5 rpca (exact 16): ", sum(conc$conc_pt15, na.rm=TRUE), "/", n, "\n")
cat("PT-1 only (exact 16):   ", sum(conc$conc_pt1,  na.rm=TRUE), "/", n, "\n")

write.xlsx(conc, file.path(output_dir, "PB90_concordance_summary.xlsx"))
cat("\n\u2713 PB90 analysis complete. Outputs in:", output_dir, "\n")