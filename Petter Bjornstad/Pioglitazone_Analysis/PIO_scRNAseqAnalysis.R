#### Pioglitazone scRNAseq Analysis - Updated Methodology
#### Uses Seurat object and celltype definitions from Sex-based analysis

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
library(ggbreak)
library(gt)
library(gtsummary)
library(fgsea)
library(msigdbr)
library(clusterProfiler)
library('org.Hs.eg.db')
library(stringr)
library(gridExtra)
library(grid)

# ============================================================================
# SETUP
# ============================================================================

setwd('C:/Users/netio/Documents/Harmonized_data/')
dir.results <- 'C:/Users/netio/Documents/UofW/Projects/pioglitazone/'

# Create output directories
dir.create(paste0(dir.results, 'AllCellTypes/'), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(dir.results, 'CellTypeSpecific/'), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(dir.results, 'GSEA/'), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(dir.results, 'VolcanoPlots/'), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(dir.results, 'Top20_DEGs/'), recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# DATA LOADING - Use Sex-based Seurat object
# ============================================================================

# Load the same Seurat object used in sex-based analysis
load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')
so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

# Apply same celltype definitions from sex-based script
so_subset$celltype1 <- case_when(
  grepl("PT-", so_subset$celltype_rpca) ~ "PT",
  grepl("TAL-", so_subset$celltype_rpca) ~ "TAL",
  grepl("EC-", so_subset$celltype_rpca) ~ "EC",
  grepl("POD", so_subset$celltype_rpca) ~ "POD",
  grepl("MAC", so_subset$celltype_rpca) ~ "MAC",
  grepl("MON", so_subset$celltype_rpca) ~ "MON",
  grepl("PC-", so_subset$celltype_rpca) ~ "PC",
  grepl("FIB", so_subset$celltype_rpca) ~ "FIB_MC_VSMC",
  grepl("DTL", so_subset$celltype_rpca) ~ "DTL",
  so_subset$celltype_rpca == "DCT" ~ "DCT",
  so_subset$celltype_rpca == "ATL" ~ "ATL",
  so_subset$celltype_rpca == "B" ~ "B",
  so_subset$celltype_rpca == "T" ~ "T"
)
so_subset$celltype1 <- as.character(so_subset$celltype1)

so_subset$KPMP_celltype2 <- as.character(so_subset$KPMP_celltype)
so_subset$celltype2 <- ifelse(
  so_subset$KPMP_celltype == "aPT" |
    so_subset$KPMP_celltype == "PT-S1/S2" |
    so_subset$KPMP_celltype == "PT-S3", "PT",
  ifelse(grepl("TAL", so_subset$KPMP_celltype), "TAL",
         ifelse(grepl("EC-", so_subset$KPMP_celltype), "EC", so_subset$KPMP_celltype2))
)

so_subset$DCT_celltype <- ifelse(
  (so_subset$KPMP_celltype == "DCT" | so_subset$KPMP_celltype == "dDCT"), "DCT", "Non-DCT"
)

# Remove CRC-55 and filter to T2D only (same as sex-based script)
so_subset <- subset(so_subset, subset = record_id != 'CRC-55')
so_subset <- subset(so_subset, subset = group == 'Type_2_Diabetes')

# ============================================================================
# MEDICATION DATA & PIOGLITAZONE GROUP ASSIGNMENT
# ============================================================================

harmonized_data <- data.table::fread("harmonized_dataset.csv")
harmonized_data <- harmonized_data %>% filter(group == 'Type 2 Diabetes')

medications <- readxl::read_xlsx("Biopsies_w_mrn_Oct3.xlsx")
medications <- medications %>% dplyr::select(mrn, ends_with('_1'), -starts_with('ever_'))
names(medications) <- str_replace(names(medications), pattern = '_1', replacement = '')

meta.data <- so_subset@meta.data

harmonized_data <- harmonized_data %>%
  filter(rh_id %in% meta.data$record_id | croc_id %in% meta.data$record_id |
           improve_id %in% meta.data$record_id | penguin_id %in% meta.data$record_id |
           rh2_id %in% meta.data$record_id) %>%
  dplyr::select(record_id, mrn, group)

medications <- medications %>% dplyr::select(mrn, tzd) %>% filter(mrn %in% harmonized_data$mrn)
harmonized_data$mrn <- as.character(harmonized_data$mrn)

final_df <- medications %>% left_join(harmonized_data, by = 'mrn')
final_df$combined_id <- paste0(final_df$mrn, '_', final_df$record_id)
final_df <- final_df %>% filter(!duplicated(combined_id))

# ============================================================================
# QUALITY CONTROL & GROUP ASSIGNMENT
# ============================================================================

scrna_small <- subset(x = so_subset, record_id %in% final_df$record_id)
scrna_small <- subset(scrna_small,
                      subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

# Assign pioglitazone groups
meta.data <- scrna_small@meta.data
new_df <- data.frame(ID = meta.data$record_id)
group_data <- final_df %>% dplyr::select(ID = record_id, group = tzd)
new_df <- new_df %>% left_join(group_data)
scrna_small$group_labels <- new_df$group

# ============================================================================
# COMPUTE POOLED OFFSET ON FULL OBJECT (before any subsetting)
# ============================================================================

counts_full <- round(GetAssayData(scrna_small, layer = 'counts'))
scrna_small$library_size <- Matrix::colSums(counts_full)
sce_full <- SingleCellExperiment(assays = list(counts = counts_full))
sce_full <- computeSumFactors(sce_full)
scrna_small$pooled_offset <- sizeFactors(sce_full)
remove(sce_full, counts_full)

# ============================================================================
# DEMOGRAPHICS TABLE
# ============================================================================

harmonized_full <- read.csv("harmonized_dataset.csv", na = '')

dat <- harmonized_full %>%
  dplyr::select(-dob) %>%
  arrange(date_of_screen) %>%
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm = TRUE))),
    .by = c(record_id, visit)
  )

record_ids_in_scrna <- unique(scrna_small$record_id)
dat_small <- dat %>% filter(record_id %in% record_ids_in_scrna)

dat_small <- dat_small %>%
  left_join(final_df %>% dplyr::select(record_id, tzd) %>% distinct(), by = 'record_id')

desc_table <- dat_small %>%
  select(age, sex, race_ethnicity, bmi, hba1c, study, tzd) %>%
  tbl_summary(
    by = tzd,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous",
      hba1c ~ "continuous",
      race_ethnicity ~ "categorical",
      study ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1, bmi ~ 1, hba1c ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      study ~ "Study"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(all_continuous() ~ "t.test")) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Pioglitazone Use**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous; n (%) for categorical")

desc_table %>%
  as_gt() %>%
  tab_options(table.font.size = 11, heading.title.font.size = 14, column_labels.font.size = 12) %>%
  gtsave(paste0(dir.results, "PIO_demographics.png"), vwidth = 1200, vheight = 800)

# ============================================================================
# UMAP VISUALIZATION
# ============================================================================

png(paste0(dir.results, 'AllCellTypes/PIO_UMAP.png'), width = 1200, height = 1200)
DimPlot(scrna_small, reduction = 'umap.harmony', group.by = 'celltype2') +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  )
dev.off()

# ============================================================================
# CELL TYPE PROPORTION ANALYSIS
# ============================================================================

cell_data <- data.frame(
  cell_type = scrna_small$celltype2,
  group = scrna_small$group_labels
) %>% filter(!is.na(group))

prop_data <- cell_data %>%
  group_by(group, cell_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(group) %>%
  mutate(total = sum(count), proportion = count / total * 100)

stat_results <- data.frame()
for (ct in unique(cell_data$cell_type)) {
  contingency <- table(cell_data$group, cell_data$cell_type == ct)
  test <- fisher.test(contingency)
  stat_results <- rbind(stat_results, data.frame(
    cell_type = ct, p_value = test$p.value
  ))
}

stat_results$p_adj <- p.adjust(stat_results$p_value, method = "BH")
stat_results$sig_label <- case_when(
  stat_results$p_adj < 0.001 ~ "***",
  stat_results$p_adj < 0.01 ~ "**",
  stat_results$p_adj < 0.05 ~ "*",
  TRUE ~ ""
)

write.csv(stat_results, paste0(dir.results, 'AllCellTypes/CellType_Proportion_FisherTests.csv'),
          row.names = FALSE)

# ============================================================================
# NEBULA ANALYSIS FUNCTION (updated for celltype2 from sex-based object)
# ============================================================================

PIO_NEBULA_Analysis <- function(so_obj, dir.results, celltype, use_hvgs = TRUE) {
  
  # Subset by cell type using celltype2 (from sex-based script)
  if (celltype == 'All') {
    so_celltype <- so_obj
  } else if (celltype %in% c('TAL', 'EC', 'PT')) {
    # celltype2 already has collapsed PT, TAL, EC
    so_celltype <- subset(so_obj, celltype2 == celltype)
  } else if (celltype == 'IC') {
    # Immune cells
    so_celltype <- subset(so_obj, celltype1 %in% c("B", "T", "MON", "MAC"))
  } else if (celltype == 'DCT') {
    so_celltype <- subset(so_obj, DCT_celltype == "DCT")
  } else if (celltype == 'POD') {
    so_celltype <- subset(so_obj, celltype1 == "POD")
  } else {
    # Try exact match on celltype2, then celltype1
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
  
  celltype2 <- str_replace_all(celltype, "/", "_")
  celltype2 <- str_replace_all(celltype2, "-", "_")
  
  cat(paste0("\n=== ", celltype2, " is running ===\n"))
  
  # Optionally filter to HVGs
  if (use_hvgs) {
    so_celltype <- FindVariableFeatures(so_celltype, selection.method = "vst", nfeatures = 2000)
    hvgs <- VariableFeatures(so_celltype)
    so_celltype <- subset(so_celltype, features = hvgs)
    prefix <- "HVG_NEBULA_"
  } else {
    prefix <- "Full_NEBULA_"
  }
  
  counts_path <- round(GetAssayData(so_celltype, layer = "counts"))
  meta_gene <- so_celltype@meta.data
  
  # Filter complete cases
  complete_idx <- complete.cases(meta_gene$group_labels)
  cat("Cells with complete data:", sum(complete_idx), "\n")
  
  meta_gene <- meta_gene[complete_idx, ]
  counts_path <- counts_path[, complete_idx]
  
  num_cells <- nrow(meta_gene)
  num_part <- length(unique(meta_gene$record_id))
  tmp_df <- meta_gene %>% dplyr::select(record_id, group_labels) %>% filter(!duplicated(record_id))
  num_yes <- tmp_df %>% filter(group_labels == 'Yes') %>% nrow()
  num_no <- tmp_df %>% filter(group_labels == 'No') %>% nrow()
  
  cat("Participants - Yes:", num_yes, "No:", num_no, "\n")
  
  # Design matrix
  pred_gene <- model.matrix(~group_labels, data = meta_gene)
  
  # Use pooled offset
  lib <- meta_gene$pooled_offset
  data_g_gene <- group_cell(count = counts_path, id = meta_gene$record_id,
                            pred = pred_gene, offset = lib)
  
  if (is.null(data_g_gene)) {
    data_g_gene <- list(count = counts_path, id = meta_gene$record_id,
                        pred = pred_gene, offset = lib)
  }
  
  # Run NEBULA
  result <- nebula(count = data_g_gene$count, id = data_g_gene$id,
                   pred = data_g_gene$pred, ncore = 1, reml = TRUE,
                   model = "NBLMM", output_re = TRUE, covariance = TRUE,
                   offset = data_g_gene$offset)
  
  full_results <- as.data.frame(result)
  full_results$num_cells <- num_cells
  full_results$num_pio_yes <- num_yes
  full_results$num_pio_no <- num_no
  
  output_file <- paste0(dir.results, prefix, celltype2, "_PIO_pooledoffset.csv")
  write.table(full_results, output_file, row.names = FALSE, quote = FALSE, sep = ',')
  
  cat(paste0(celltype2, " is done. Saved to: ", output_file, "\n"))
  return(invisible(full_results))
}

# ============================================================================
# RUN NEBULA FOR ALL CELL TYPES (HVG version)
# ============================================================================

celltypes_vec <- c('All', 'PT', 'TAL', 'EC', 'POD', 'DCT', 'IC')

for (i in 1:length(celltypes_vec)) {
  PIO_NEBULA_Analysis(so_obj = scrna_small, dir.results = paste0(dir.results, 'CellTypeSpecific/'),
                      celltype = celltypes_vec[i], use_hvgs = TRUE)
  cat(paste0("Completed: ", celltypes_vec[i], "\n"))
}

# ============================================================================
# RUN NEBULA FOR ALL CELL TYPES (Full gene version)
# ============================================================================

for (i in 1:length(celltypes_vec)) {
  PIO_NEBULA_Analysis(so_obj = scrna_small, dir.results = paste0(dir.results, 'CellTypeSpecific/'),
                      celltype = celltypes_vec[i], use_hvgs = FALSE)
  cat(paste0("Completed (full): ", celltypes_vec[i], "\n"))
}

# ============================================================================
# VOLCANO PLOTS (with axis break logic)
# ============================================================================

for (i in 1:length(celltypes_vec)) {
  
  ct <- celltypes_vec[i]
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  
  file_path <- paste0(dir.results, 'CellTypeSpecific/Full_NEBULA_', ct_clean, '_PIO_pooledoffset.csv')
  if (!file.exists(file_path)) next
  
  sig_markers <- data.table::fread(file_path)
  sig_markers <- sig_markers %>% dplyr::select(
    Gene = summary.gene,
    LogFC = summary.logFC_group_labelsYes,
    Pvalue = summary.p_group_labelsYes
  )
  
  tmp_df <- sig_markers %>%
    mutate(
      p_adj = p.adjust(Pvalue, method = 'BH'),
      diffexp = case_when(
        Pvalue < 0.05 & LogFC > 0 ~ 'Up',
        Pvalue < 0.05 & LogFC < 0 ~ 'Down',
        TRUE ~ 'No'
      )
    ) %>%
    arrange(Pvalue) %>%
    mutate(label = ifelse(row_number() <= 10, Gene, NA)) %>%
    filter(abs(LogFC) < 10)
  
  # Axis break detection
  logfc_q95 <- quantile(abs(tmp_df$LogFC), 0.95, na.rm = TRUE)
  logfc_max <- max(abs(tmp_df$LogFC), na.rm = TRUE)
  needs_break <- logfc_max > (logfc_q95 * 2)
  
  use_break <- FALSE
  if (needs_break) {
    sorted_logfc <- sort(abs(tmp_df$LogFC))
    gaps <- diff(sorted_logfc)
    gap_threshold <- quantile(abs(tmp_df$LogFC), 0.90, na.rm = TRUE)
    outer_data_idx <- which(sorted_logfc > gap_threshold)
    if (length(outer_data_idx) > 1) {
      outer_gaps <- gaps[outer_data_idx[-length(outer_data_idx)]]
      max_gap_idx <- which.max(outer_gaps)
      break_start <- sorted_logfc[outer_data_idx[max_gap_idx]] + 0.1
      break_end <- sorted_logfc[outer_data_idx[max_gap_idx + 1]] - 0.1
      if (break_end - break_start > 0.5) use_break <- TRUE
    }
  }
  
  if (length(unique(tmp_df$diffexp)) > 1) {
    tmp_graph <- ggplot(tmp_df, aes(x = LogFC, y = -log10(Pvalue), col = diffexp, label = label)) +
      geom_point() +
      geom_text(size = 2, vjust = 2, color = 'black') +
      scale_color_manual(values = c('orange', 'grey', 'purple'),
                         labels = c('Downregulated', 'Not significant', 'Upregulated')) +
      geom_hline(yintercept = -log10(0.05), col = 'blue', linetype = 'dashed') +
      geom_vline(xintercept = 0, col = 'black', linetype = 'dashed') +
      theme_classic() +
      labs(x = 'LogFC', y = '-log10 P-value', col = 'Differential Expression',
           title = paste0('Pioglitazone Effects in ', ct, ' Cells'))
  } else {
    tmp_graph <- ggplot(tmp_df, aes(x = LogFC, y = -log10(Pvalue), col = diffexp, label = label)) +
      geom_point() +
      geom_text(size = 2, vjust = 2, color = 'black') +
      scale_color_manual(values = c('grey'), labels = c('Not significant')) +
      geom_hline(yintercept = -log10(0.05), col = 'blue', linetype = 'dashed') +
      geom_vline(xintercept = 0, col = 'black', linetype = 'dashed') +
      theme_classic() +
      labs(x = 'LogFC', y = '-log10 P-value', col = 'Differential Expression',
           title = paste0('Pioglitazone Effects in ', ct, ' Cells'))
  }
  
  if (use_break) {
    tmp_graph <- tmp_graph + scale_x_break(breaks = c(break_start, break_end), scales = 0.3)
  }
  
  pdf(paste0(dir.results, 'VolcanoPlots/Volcano_', ct_clean, '_PIO.pdf'))
  print(tmp_graph)
  dev.off()
  
  png(paste0(dir.results, 'VolcanoPlots/Volcano_', ct_clean, '_PIO.png'),
      width = 3000, height = 3000, res = 300)
  print(tmp_graph)
  dev.off()
  
  cat(paste0('Volcano plot done for ', ct, if (use_break) ' (with axis break)' else '', '\n'))
}

# ============================================================================
# COMBINED 4-PANEL VOLCANO PLOT
# ============================================================================

plot_celltypes <- c('All', 'PT', 'TAL', 'EC')
plot_list <- list()
all_logfc <- c(); all_pval <- c()

# First pass: get common axis limits
for (i in 1:length(plot_celltypes)) {
  ct_clean <- str_replace_all(str_replace_all(plot_celltypes[i], "/", "_"), "-", "_")
  fp <- paste0(dir.results, 'CellTypeSpecific/Full_NEBULA_', ct_clean, '_PIO_pooledoffset.csv')
  if (!file.exists(fp)) next
  d <- data.table::fread(fp) %>%
    dplyr::select(LogFC = summary.logFC_group_labelsYes, Pvalue = summary.p_group_labelsYes) %>%
    filter(abs(LogFC) < 10)
  all_logfc <- c(all_logfc, d$LogFC)
  all_pval <- c(all_pval, -log10(d$Pvalue))
}

x_limit <- max(abs(all_logfc), na.rm = TRUE) * 1.05
y_limit <- max(all_pval, na.rm = TRUE) * 1.05

# Second pass: create plots
for (i in 1:length(plot_celltypes)) {
  ct <- plot_celltypes[i]
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  fp <- paste0(dir.results, 'CellTypeSpecific/Full_NEBULA_', ct_clean, '_PIO_pooledoffset.csv')
  if (!file.exists(fp)) next
  
  d <- data.table::fread(fp) %>%
    dplyr::select(Gene = summary.gene, LogFC = summary.logFC_group_labelsYes,
                  Pvalue = summary.p_group_labelsYes) %>%
    mutate(diffexp = case_when(
      Pvalue < 0.05 & LogFC > 0 ~ 'Up',
      Pvalue < 0.05 & LogFC < 0 ~ 'Down',
      TRUE ~ 'No'
    )) %>%
    arrange(Pvalue) %>%
    mutate(label = ifelse(row_number() <= 10, Gene, NA)) %>%
    filter(abs(LogFC) < 10)
  
  color_vals <- c('orange', 'grey', 'purple')[c('Down', 'No', 'Up') %in% unique(d$diffexp)]
  color_labs <- c('Down in Pio', 'Not significant', 'Up in Pio')[c('Down', 'No', 'Up') %in% unique(d$diffexp)]
  
  plot_list[[i]] <- ggplot(d, aes(x = LogFC, y = -log10(Pvalue), col = diffexp, label = label)) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_text(size = 3, vjust = 2, color = 'black') +
    scale_color_manual(values = setNames(color_vals, unique(d$diffexp)[order(unique(d$diffexp))]),
                       labels = color_labs) +
    geom_hline(yintercept = -log10(0.05), col = 'blue', linetype = 'dashed', linewidth = 0.5) +
    geom_vline(xintercept = 0, col = 'black', linetype = 'dashed', linewidth = 0.5) +
    coord_cartesian(xlim = c(-x_limit, x_limit), ylim = c(0, y_limit)) +
    theme_classic() +
    labs(x = 'LogFC (Pio Yes vs No)', y = '-log10 P-value', col = 'Differential Expression',
         title = paste0(ct, ' Cells'), tag = LETTERS[i]) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          plot.tag = element_text(size = 16, face = "bold"),
          aspect.ratio = 1)
}

combined_volcano <- (plot_list[[1]] | plot_list[[2]]) /
  (plot_list[[3]] | plot_list[[4]]) +
  plot_layout(guides = 'collect')

ggsave(paste0(dir.results, 'VolcanoPlots/Combined_4Panel_PIO.pdf'),
       combined_volcano, width = 14, height = 12)
ggsave(paste0(dir.results, 'VolcanoPlots/Combined_4Panel_PIO.png'),
       combined_volcano, width = 14, height = 12, dpi = 300)

# ============================================================================
# TOP 20 DEG TABLES WITH FDR
# ============================================================================

for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  fp <- paste0(dir.results, 'CellTypeSpecific/Full_NEBULA_', ct_clean, '_PIO_pooledoffset.csv')
  if (!file.exists(fp)) next
  
  data <- read.csv(fp) %>% filter(!is.na(summary.p_group_labelsYes))
  
  cat(paste("\n---", ct, "---\n"))
  cat("Total genes:", nrow(data), "\n")
  cat("P < 0.05:", sum(data$summary.p_group_labelsYes < 0.05), "\n")
  
  data$adjusted_pvalue <- p.adjust(data$summary.p_group_labelsYes, method = "BH")
  cat("FDR < 0.05:", sum(data$adjusted_pvalue < 0.05, na.rm = TRUE), "\n")
  
  top20 <- data %>%
    arrange(summary.p_group_labelsYes) %>%
    head(20) %>%
    transmute(
      gene = summary.gene,
      pvalue = summary.p_group_labelsYes,
      adjusted_pvalue = adjusted_pvalue,
      logFC = summary.logFC_group_labelsYes,
      cell_number = num_cells
    )
  
  write.csv(top20, paste0(dir.results, 'Top20_DEGs/Top20_', ct_clean, '_PIO.csv'), row.names = FALSE)
}

# ============================================================================
# GSEA (Hallmark + GO:BP, GO:CC, GO:MF)
# ============================================================================

hallmark_list <- split(msigdbr(species = "Homo sapiens", category = "H")$gene_symbol,
                       msigdbr(species = "Homo sapiens", category = "H")$gs_name)
go_bp_list <- split(msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")$gene_symbol,
                    msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")$gs_name)
go_cc_list <- split(msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")$gene_symbol,
                    msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")$gs_name)
go_mf_list <- split(msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")$gene_symbol,
                    msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")$gs_name)

geneset_types <- list(Hallmark = hallmark_list, GO_BP = go_bp_list,
                      GO_CC = go_cc_list, GO_MF = go_mf_list)

all_gsea_results <- list()

for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  fp <- paste0(dir.results, 'CellTypeSpecific/Full_NEBULA_', ct_clean, '_PIO_pooledoffset.csv')
  if (!file.exists(fp)) next
  
  cat("\n=== Processing GSEA for", ct, "===\n")
  
  de_results <- read.csv(fp)
  ranked_genes <- de_results %>%
    dplyr::select(logFC = summary.logFC_group_labelsYes, gene_id = summary.gene) %>%
    filter(!is.na(logFC)) %>%
    arrange(desc(logFC)) %>%
    pull(logFC, name = gene_id)
  
  for (gs_name in names(geneset_types)) {
    cat("  Running fgsea for", gs_name, "...\n")
    
    gsea_res <- fgsea(pathways = geneset_types[[gs_name]], stats = ranked_genes,
                      minSize = 15, maxSize = 500, nPermSimple = 10000)
    
    gsea_res$celltype <- ct
    gsea_res$geneset_type <- gs_name
    
    result_key <- paste0(ct, "_", gs_name)
    all_gsea_results[[result_key]] <- gsea_res
    
    top_paths <- gsea_res %>%
      filter(pval < 0.05) %>%
      arrange(pval) %>%
      slice_head(n = 20) %>%
      mutate(pathway_clean = gsub("HALLMARK_|GOBP_|GOCC_|GOMF_", "", pathway),
             pathway_clean = gsub("_", " ", pathway_clean),
             pathway_clean = tools::toTitleCase(tolower(pathway_clean)))
    
    if (nrow(top_paths) > 0) {
      p <- ggplot(top_paths, aes(x = NES, y = reorder(pathway_clean, NES))) +
        geom_point(aes(size = -log10(pval), color = NES)) +
        scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
        scale_size_continuous(range = c(2, 10)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        theme_bw() +
        theme(axis.text.y = element_text(size = 9),
              plot.title = element_text(face = "bold", hjust = 0.5)) +
        labs(x = "NES", y = "", title = paste0(ct, " - ", gs_name, " (Top 20)"))
      
      ggsave(paste0(dir.results, 'GSEA/', gs_name, '_gsea_PIO_', ct_clean, '.pdf'),
             plot = p, width = 10, height = 8)
      ggsave(paste0(dir.results, 'GSEA/', gs_name, '_gsea_PIO_', ct_clean, '.png'),
             plot = p, width = 10, height = 8, dpi = 300)
    }
  }
}

# Save combined GSEA results
combined_gsea <- bind_rows(all_gsea_results)
combined_gsea_save <- combined_gsea %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";")))
write.csv(combined_gsea_save, paste0(dir.results, 'GSEA/all_celltypes_gsea_PIO_results.csv'),
          row.names = FALSE)

for (gs_name in names(geneset_types)) {
  subset_res <- combined_gsea %>%
    filter(geneset_type == gs_name) %>%
    mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";")))
  write.csv(subset_res, paste0(dir.results, 'GSEA/gsea_results_PIO_', gs_name, '.csv'),
            row.names = FALSE)
}

# Summary table
summary_table <- combined_gsea %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";"))) %>%
  filter(pval < 0.05) %>%
  group_by(celltype, geneset_type) %>%
  summarise(n_significant = n(), .groups = 'drop') %>%
  pivot_wider(names_from = geneset_type, values_from = n_significant, values_fill = 0)

print(summary_table)
write.csv(summary_table, paste0(dir.results, 'GSEA/summary_significant_pathways_PIO.csv'),
          row.names = FALSE)

cat("\n=== Pioglitazone analysis complete! ===\n")
