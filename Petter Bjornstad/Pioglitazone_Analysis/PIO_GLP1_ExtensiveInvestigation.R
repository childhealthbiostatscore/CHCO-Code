# ============================================================================
# GLP-1 x PIOGLITAZONE EXTENDED ANALYSIS
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
# DIRECTORIES
# ============================================================================

setwd('C:/Users/netio/Documents/Harmonized_data/')

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
# LOAD SEURAT OBJECT
# ============================================================================

load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')

# Re-set dirs immediately after load() since RData may overwrite them
dir.sensitivity <- 'C:/Users/netio/Documents/UofW/Projects/pioglitazone/Sensitivity/'
dir.interaction <- 'C:/Users/netio/Documents/UofW/Projects/pioglitazone/Sensitivity/Interaction_PIO_GLP1/'
dir.root        <- 'C:/Users/netio/Documents/UofW/Projects/pioglitazone/GLP1_Extended/'

so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

# ---- Cell type definitions (exact copy from sensitivity script) -------------
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
         ifelse(grepl("EC-", so_subset$KPMP_celltype), "EC",
                so_subset$KPMP_celltype2))
)

so_subset$DCT_celltype <- ifelse(
  so_subset$KPMP_celltype %in% c("DCT", "dDCT"), "DCT", "Non-DCT"
)

# ---- Filter ----------------------------------------------------------------
so_subset <- subset(so_subset, subset = record_id != 'CRC-55')
so_subset <- subset(so_subset, subset = group == 'Type_2_Diabetes')

# ============================================================================
# MEDICATION & COVARIATE DATA
# (exact same approach as sensitivity script)
# ============================================================================

harmonized_data <- read.csv(
  "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv",
  na = '')

harmonized_data <- harmonized_data %>% filter(group == 'Type 2 Diabetes')

harmonized_data <- harmonized_data %>%
  dplyr::select(-dob) %>%
  arrange(date_of_screen) %>%
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric),         ~ ifelse(all(is.na(.x)), NA_real_,      mean(na.omit(.x), na.rm = TRUE))),
    .by = c(record_id, visit)
  )

medications <- readxl::read_xlsx("Biopsies_w_mrn_Oct3.xlsx")
medications <- medications %>% dplyr::select(mrn, ends_with('_1'), -starts_with('ever_'))
names(medications) <- str_replace(names(medications), pattern = '_1', replacement = '')

meta.data <- so_subset@meta.data

harmonized_data <- harmonized_data %>%
  filter(rh_id %in% meta.data$record_id | croc_id %in% meta.data$record_id |
           improve_id %in% meta.data$record_id | penguin_id %in% meta.data$record_id |
           rh2_id %in% meta.data$record_id) %>%
  dplyr::select(record_id, mrn, group, bmi, sex, age, race_ethnicity,
                hba1c, eGFR_CKD_epi, acr_u,
                epic_sglti2_1, epic_insulin_1, epic_raasi_1,
                epic_statin_1, epic_fibrate_1)

medications <- medications %>%
  dplyr::select(mrn, tzd, epic_mfm_1 = mfm, epic_glp1ra_1 = glp1ra) %>%
  filter(mrn %in% harmonized_data$mrn)

harmonized_data$mrn <- as.character(harmonized_data$mrn)

final_df <- medications %>% left_join(harmonized_data, by = 'mrn')
final_df$combined_id <- paste0(final_df$mrn, '_', final_df$record_id)
final_df <- final_df %>% filter(!duplicated(combined_id))

# ============================================================================
# QC & APPLY GROUP LABELS TO SEURAT OBJECT
# ============================================================================

scrna_small <- subset(so_subset, record_id %in% final_df$record_id)
scrna_small <- subset(scrna_small,
                      subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
rm(so_subset)

# Drop any conflicting covariate columns before join (prevents .x/.y clobbering)
meta.data <- scrna_small@meta.data %>%
  dplyr::select(-any_of(c("bmi", "sex", "epic_mfm_1", "epic_glp1ra_1", "group_labels",
                          "age", "race_ethnicity", "hba1c", "eGFR_CKD_epi", "acr_u",
                          "epic_sglti2_1", "epic_insulin_1", "epic_raasi_1",
                          "epic_statin_1", "epic_fibrate_1")))

covar_lookup <- final_df %>%
  dplyr::select(record_id, tzd, bmi, sex, epic_mfm_1, epic_glp1ra_1,
                age, race_ethnicity, hba1c, eGFR_CKD_epi, acr_u,
                epic_sglti2_1, epic_insulin_1, epic_raasi_1,
                epic_statin_1, epic_fibrate_1) %>%
  distinct(record_id, .keep_all = TRUE)

meta_joined <- meta.data %>%
  left_join(covar_lookup, by = 'record_id')

scrna_small$group_labels   <- meta_joined$tzd
scrna_small$bmi            <- meta_joined$bmi
scrna_small$sex            <- meta_joined$sex
scrna_small$epic_mfm_1     <- meta_joined$epic_mfm_1
scrna_small$epic_glp1ra_1  <- meta_joined$epic_glp1ra_1
scrna_small$age            <- meta_joined$age
scrna_small$race_ethnicity <- meta_joined$race_ethnicity
scrna_small$hba1c          <- meta_joined$hba1c
scrna_small$eGFR_CKD_epi   <- meta_joined$eGFR_CKD_epi
scrna_small$acr_u          <- meta_joined$acr_u
scrna_small$epic_sglti2_1  <- meta_joined$epic_sglti2_1
scrna_small$epic_insulin_1 <- meta_joined$epic_insulin_1
scrna_small$epic_raasi_1   <- meta_joined$epic_raasi_1
scrna_small$epic_statin_1  <- meta_joined$epic_statin_1
scrna_small$epic_fibrate_1 <- meta_joined$epic_fibrate_1

cat("group_labels (PIO):\n");    print(table(scrna_small$group_labels,  useNA = "ifany"))
cat("epic_glp1ra_1 (GLP-1RA):\n"); print(table(scrna_small$epic_glp1ra_1, useNA = "ifany"))

# ============================================================================
# POOLED OFFSET
# ============================================================================

counts_full <- round(GetAssayData(scrna_small, layer = 'counts'))
scrna_small$library_size <- Matrix::colSums(counts_full)
sce_full <- SingleCellExperiment(assays = list(counts = counts_full))
sce_full <- computeSumFactors(sce_full)
scrna_small$pooled_offset <- sizeFactors(sce_full)
rm(sce_full, counts_full)

# ============================================================================
# SHARED CONSTANTS
# ============================================================================

celltypes_vec   <- c('All', 'PT', 'TAL', 'EC', 'POD', 'DCT', 'IC')
pt_subtypes_vec <- c('aPT', 'PT-S1/S2', 'PT-S3')
GSEA_PARAMS     <- list(minSize = 15, maxSize = 500, nPermSimple = 10000)
P_THRESH        <- 0.05

GROUP_COLORS <- c("PIO-GLP1-" = "#999999", "PIO+GLP1-" = "#E69F00",
                  "PIO-GLP1+" = "#56B4E9", "PIO+GLP1+" = "#009E73")

hallmark_list <- split(msigdbr(species = "Homo sapiens", category = "H")$gene_symbol,
                       msigdbr(species = "Homo sapiens", category = "H")$gs_name)
go_bp_list    <- split(msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")$gene_symbol,
                       msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")$gs_name)
geneset_types <- list(Hallmark = hallmark_list, GO_BP = go_bp_list)

# ============================================================================
# FOUR-GROUP LABEL
# ============================================================================

scrna_small$four_group <- case_when(
  scrna_small$group_labels  == "Yes" & scrna_small$epic_glp1ra_1 == "Yes" ~ "PIO+GLP1+",
  scrna_small$group_labels  == "Yes" & scrna_small$epic_glp1ra_1 == "No"  ~ "PIO+GLP1-",
  scrna_small$group_labels  == "No"  & scrna_small$epic_glp1ra_1 == "Yes" ~ "PIO-GLP1+",
  scrna_small$group_labels  == "No"  & scrna_small$epic_glp1ra_1 == "No"  ~ "PIO-GLP1-",
  TRUE ~ NA_character_
)
scrna_small$four_group <- factor(scrna_small$four_group,
                                 levels = c("PIO-GLP1-", "PIO+GLP1-", "PIO-GLP1+", "PIO+GLP1+"))

cat("Four-group distribution:\n")
print(table(scrna_small$four_group, useNA = "ifany"))

# ============================================================================
# SECTION 1 — FOUR-GROUP DEMOGRAPHICS
# ============================================================================

cat("\n=== SECTION 1: Four-Group Demographics ===\n")

patient_meta <- scrna_small@meta.data %>%
  as.data.frame() %>%
  mutate(four_group = as.character(four_group)) %>%
  filter(!duplicated(record_id)) %>%
  filter(!is.na(four_group)) %>%
  mutate(four_group = factor(four_group, levels = c("PIO-GLP1-", "PIO+GLP1-", "PIO-GLP1+", "PIO+GLP1+")))

group_counts <- patient_meta %>%
  count(four_group, name = "n_patients") %>%
  mutate(pct = round(100 * n_patients / sum(n_patients), 1))
cat("\nPatient counts:\n"); print(group_counts)
write.csv(group_counts, paste0(dir.root, 'Demographics/PatientCounts_FourGroup.csv'), row.names = FALSE)

# ---- Full gtsummary demographics table -------------------------------------

demo_vars <- c("age", "sex", "race_ethnicity", "bmi", "hba1c", "eGFR_CKD_epi", "acr_u",
               "epic_glp1ra_1", "epic_mfm_1", "epic_sglti2_1",
               "epic_insulin_1", "epic_raasi_1", "epic_statin_1", "epic_fibrate_1")

demo_vars_present <- intersect(demo_vars, colnames(patient_meta))
missing_vars      <- setdiff(demo_vars, colnames(patient_meta))
if (length(missing_vars) > 0)
  cat("Skipping missing demo vars:", paste(missing_vars, collapse = ", "), "\n")

cont_candidates <- c("age", "bmi", "hba1c", "eGFR_CKD_epi", "acr_u")
cont_present    <- intersect(cont_candidates, demo_vars_present)
cat_present     <- setdiff(demo_vars_present, cont_present)

type_list  <- c(setNames(rep(list("continuous"),  length(cont_present)),  cont_present),
                setNames(rep(list("categorical"), length(cat_present)),   cat_present))

label_list <- list(
  age            ~ "Age, years",
  sex            ~ "Sex",
  race_ethnicity ~ "Race/Ethnicity",
  bmi            ~ "BMI, kg/m²",
  hba1c          ~ "HbA1c, %",
  eGFR_CKD_epi   ~ "eGFR, mL/min/1.73m²",
  acr_u          ~ "UACR, mg/g",
  epic_glp1ra_1  ~ "GLP-1RA",
  epic_mfm_1     ~ "Metformin",
  epic_sglti2_1  ~ "SGLT2i",
  epic_insulin_1 ~ "Insulin",
  epic_raasi_1   ~ "ACEi/ARB",
  epic_statin_1  ~ "Statin",
  epic_fibrate_1 ~ "Fenofibrate"
)
label_list <- label_list[names(label_list) %in% demo_vars_present]

demo_table <- patient_meta %>%
  dplyr::select(four_group, all_of(demo_vars_present)) %>%
  tbl_summary(
    by        = four_group,
    type      = type_list,
    statistic = list(all_continuous()  ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits    = list(all_continuous()  ~ 1,
                     all_categorical() ~ c(0, 1)),
    label        = label_list,
    missing_text = "Missing"
  ) %>%
  add_p(test = list(all_continuous()  ~ "kruskal.test",
                    all_categorical() ~ "chisq.test")) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**PIO x GLP-1RA Group**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) or n (%); Kruskal-Wallis / χ²")

demo_table %>%
  as_gt() %>%
  tab_options(table.font.size = 11, heading.title.font.size = 14,
              column_labels.font.size = 12) %>%
  gtsave(paste0(dir.root, 'Demographics/FourGroup_Demographics.png'), vwidth = 1400, vheight = 1000)
demo_table %>% as_gt() %>%
  gtsave(paste0(dir.root, 'Demographics/FourGroup_Demographics.html'))

# ---- Cell counts bar plot --------------------------------------------------
cell_counts <- scrna_small@meta.data %>%
  filter(!is.na(four_group)) %>%
  count(four_group, celltype2, name = "n_cells")
write.csv(cell_counts, paste0(dir.root, 'Demographics/CellCounts_FourGroup_byCelltype.csv'), row.names = FALSE)

p_cells <- ggplot(cell_counts, aes(x = celltype2, y = n_cells, fill = four_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = GROUP_COLORS) +
  theme_bw() +
  labs(x = "Cell Type", y = "Number of Cells", fill = "Group",
       title = "Cell counts by four-group (PIO x GLP-1RA)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title  = element_text(face = "bold"))
ggsave(paste0(dir.root, 'Demographics/CellCounts_FourGroup_barplot.pdf'), p_cells, width = 10, height = 5)
ggsave(paste0(dir.root, 'Demographics/CellCounts_FourGroup_barplot.png'), p_cells, width = 10, height = 5, dpi = 200)

# ---- UMAP ------------------------------------------------------------------
if ("UMAP_1" %in% colnames(scrna_small@meta.data)) {
  umap_df <- scrna_small@meta.data %>%
    filter(!is.na(four_group)) %>%
    dplyr::select(UMAP_1, UMAP_2, four_group)
  p_umap <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = four_group)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_manual(values = GROUP_COLORS,
                       guide  = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_classic() +
    labs(color = "Group", title = "UMAP — PIO x GLP-1RA Groups") +
    theme(plot.title = element_text(face = "bold"))
  ggsave(paste0(dir.root, 'Demographics/UMAP_FourGroup.pdf'), p_umap, width = 8, height = 6)
  ggsave(paste0(dir.root, 'Demographics/UMAP_FourGroup.png'), p_umap, width = 8, height = 6, dpi = 200)
}

cat("Section 1 complete.\n")

# ============================================================================
# SECTION 2 — GLP-1RA MAIN EFFECT (NEBULA)
# ============================================================================

cat("\n=== SECTION 2: GLP-1RA Main Effect (NEBULA) ===\n")

GLP1_COVARIATES <- list(
  base      = NULL,
  plus_pio  = "group_labels",
  plus_bmi  = "bmi",
  plus_sex  = "sex",
  plus_mfm  = "epic_mfm_1"
)

run_nebula_glp1_main <- function(so_obj, dir.out, celltype,
                                 covariate_name = NULL, model_label = "base") {
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
    else { cat("Cell type", celltype, "not found. Skipping.\n"); return(invisible(NULL)) }
  }
  
  DefaultAssay(so_c) <- "RNA"
  ct_clean <- str_replace_all(str_replace_all(celltype, "/", "_"), "-", "_")
  cat(paste0("\n=== GLP1 main | ", ct_clean, " | model: ", model_label, " ===\n"))
  
  counts_mat <- round(GetAssayData(so_c, layer = "counts"))
  meta_g     <- so_c@meta.data
  
  required_cols <- c("epic_glp1ra_1", if (!is.null(covariate_name)) covariate_name)
  complete_idx  <- complete.cases(meta_g[, required_cols, drop = FALSE])
  meta_g        <- meta_g[complete_idx, ]
  counts_mat    <- counts_mat[, complete_idx]
  
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

for (ct in celltypes_vec) {
  for (m in names(GLP1_COVARIATES)) {
    run_nebula_glp1_main(so_obj = scrna_small, dir.out = dir.root,
                         celltype = ct, covariate_name = GLP1_COVARIATES[[m]],
                         model_label = m)
  }
}

# Helper: load GLP-1 effect
load_glp1_effect <- function(dir.root, ct_clean, model_label = "base") {
  fp <- paste0(dir.root, 'GLP1_Main/NEBULA/GLP1_NEBULA_', ct_clean, '_', model_label, '.csv')
  if (!file.exists(fp)) { cat("  [load_glp1] Not found:", fp, "\n"); return(NULL) }
  df <- read.csv(fp)
  
  logfc_col <- grep("logFC.*glp1ra.*Yes|logFC.*epic_glp1ra_1Yes",
                    names(df), value = TRUE, ignore.case = TRUE)[1]
  p_col     <- grep("^summary\\.p_.*glp1ra.*Yes|^summary\\.p_.*epic_glp1ra_1Yes",
                    names(df), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(logfc_col)) {
    logfc_col <- grep("logFC", names(df), value = TRUE)[2]
    p_col     <- grep("^summary\\.p_", names(df), value = TRUE)[2]
  }
  gene_col <- if ("summary.gene" %in% names(df)) "summary.gene" else names(df)[1]
  
  df %>% transmute(gene   = .data[[gene_col]],
                   logFC  = .data[[logfc_col]],
                   pvalue = .data[[p_col]],
                   fdr    = p.adjust(pvalue, method = "BH"))
}

# ---- Volcano plots ---------------------------------------------------------
for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  vol_list <- list()
  for (m in names(GLP1_COVARIATES)) {
    df <- load_glp1_effect(dir.root, ct_clean, m)
    if (is.null(df)) next
    df <- df %>%
      mutate(diffexp = case_when(pvalue < P_THRESH & logFC > 0 ~ "Up",
                                 pvalue < P_THRESH & logFC < 0 ~ "Down",
                                 TRUE ~ "NS")) %>%
      arrange(pvalue) %>%
      mutate(label = ifelse(row_number() <= 8, gene, NA)) %>%
      filter(abs(logFC) < 10)
    vol_list[[m]] <- ggplot(df, aes(x = logFC, y = -log10(pvalue), color = diffexp, label = label)) +
      geom_point(size = 1.2, alpha = 0.7) +
      ggrepel::geom_text_repel(size = 2.2, color = "black", max.overlaps = 12, na.rm = TRUE) +
      scale_color_manual(values = c(Down = "steelblue", NS = "grey70", Up = "firebrick")) +
      geom_hline(yintercept = -log10(P_THRESH), color = "blue", linetype = "dashed", linewidth = 0.4) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
      theme_classic() +
      labs(x = "LogFC", y = "-log10(p)", title = paste0(ct, ": GLP-1 ", m)) +
      theme(plot.title = element_text(face = "bold", size = 10), legend.position = "bottom", aspect.ratio = 1)
  }
  if (length(vol_list) > 1) {
    combined_vol <- wrap_plots(vol_list, ncol = 3) +
      plot_layout(guides = "collect") & theme(legend.position = "bottom")
    ggsave(paste0(dir.root, 'GLP1_Main/VolcanoPlots/GLP1_Volcano_', ct_clean, '.pdf'), combined_vol, width = 18, height = 6)
    ggsave(paste0(dir.root, 'GLP1_Main/VolcanoPlots/GLP1_Volcano_', ct_clean, '.png'), combined_vol, width = 18, height = 6, dpi = 200)
  }
}

# ---- GSEA ------------------------------------------------------------------
all_glp1_gsea <- list()
for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  for (m in names(GLP1_COVARIATES)) {
    df <- load_glp1_effect(dir.root, ct_clean, m)
    if (is.null(df)) next
    ranked <- df %>% filter(!is.na(logFC)) %>% arrange(desc(logFC)) %>%
      { setNames(.$logFC, .$gene) }
    for (gs_name in names(geneset_types)) {
      key <- paste(ct, m, gs_name, sep = "__")
      cat("  GLP1 GSEA:", key, "\n")
      set.seed(42)
      res <- fgsea(pathways = geneset_types[[gs_name]], stats = ranked,
                   minSize = GSEA_PARAMS$minSize, maxSize = GSEA_PARAMS$maxSize,
                   nPermSimple = GSEA_PARAMS$nPermSimple)
      res$celltype      <- ct; res$model <- m; res$geneset_type <- gs_name
      res$pathway_clean <- gsub("HALLMARK_|GOBP_", "", res$pathway) %>%
        gsub("_", " ", .) %>% tools::toTitleCase()
      all_glp1_gsea[[key]] <- res
    }
  }
}

combined_glp1_gsea <- bind_rows(all_glp1_gsea) %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.csv(combined_glp1_gsea, paste0(dir.root, 'GLP1_Main/GSEA/Raw/GLP1_all_gsea.csv'), row.names = FALSE)
cat("Section 2 complete.\n")

# ============================================================================
# SECTION 3 — PT SUBTYPE DEEP-DIVE
# ============================================================================

cat("\n=== SECTION 3: PT Subtype Analysis ===\n")

scrna_small$pt_subtype <- case_when(
  scrna_small$KPMP_celltype == "aPT"      ~ "aPT",
  scrna_small$KPMP_celltype == "PT-S1/S2" ~ "PT-S1/S2",
  scrna_small$KPMP_celltype == "PT-S3"    ~ "PT-S3",
  TRUE ~ NA_character_
)

observed_pt <- scrna_small@meta.data %>%
  filter(!is.na(pt_subtype)) %>%
  count(pt_subtype) %>%
  filter(n > 50)
cat("PT subtypes with >50 cells:\n"); print(observed_pt)

pt_models <- list(
  pio_base    = list(exposure = "group_labels",  covariate = NULL),
  pio_glp1adj = list(exposure = "group_labels",  covariate = "epic_glp1ra_1"),
  glp1_base   = list(exposure = "epic_glp1ra_1", covariate = NULL),
  glp1_pioadj = list(exposure = "epic_glp1ra_1", covariate = "group_labels")
)

run_nebula_pt_subtype <- function(so_obj, dir.out, pt_sub, exposure,
                                  covariate_name = NULL, model_label) {
  so_c <- subset(so_obj, pt_subtype == pt_sub)
  if (ncol(so_c) < 50) { cat("Skipping", pt_sub, "— <50 cells.\n"); return(invisible(NULL)) }
  DefaultAssay(so_c) <- "RNA"
  sub_clean  <- str_replace_all(str_replace_all(pt_sub, "/", "_"), "-", "_")
  counts_mat <- round(GetAssayData(so_c, layer = "counts"))
  meta_g     <- so_c@meta.data
  required_cols <- c(exposure, if (!is.null(covariate_name)) covariate_name)
  complete_idx  <- complete.cases(meta_g[, required_cols, drop = FALSE])
  meta_g     <- meta_g[complete_idx, ]; counts_mat <- counts_mat[, complete_idx]
  formula_str <- if (is.null(covariate_name)) paste0("~", exposure) else
    paste0("~", exposure, " + ", covariate_name)
  pred_gene <- model.matrix(as.formula(formula_str), data = meta_g)
  lib    <- meta_g$pooled_offset
  data_g <- group_cell(count = counts_mat, id = meta_g$record_id, pred = pred_gene, offset = lib)
  if (is.null(data_g)) data_g <- list(count = counts_mat, id = meta_g$record_id, pred = pred_gene, offset = lib)
  result <- nebula(count = data_g$count, id = data_g$id, pred = data_g$pred,
                   ncore = 1, reml = TRUE, model = "NBLMM", output_re = TRUE,
                   covariance = TRUE, offset = data_g$offset)
  out_df <- as.data.frame(result)
  out_df$pt_subtype  <- pt_sub; out_df$exposure <- exposure
  out_df$model_label <- model_label; out_df$num_cells <- nrow(meta_g)
  out_file <- paste0(dir.out, 'PT_Subtypes/NEBULA/PT_NEBULA_', sub_clean, '_', model_label, '.csv')
  write.csv(out_df, out_file, row.names = FALSE)
  cat("Saved:", out_file, "\n")
  return(invisible(out_df))
}

for (pt_sub in observed_pt$pt_subtype) {
  for (m in names(pt_models)) {
    run_nebula_pt_subtype(so_obj = scrna_small, dir.out = dir.root, pt_sub = pt_sub,
                          exposure = pt_models[[m]]$exposure,
                          covariate_name = pt_models[[m]]$covariate, model_label = m)
  }
}
cat("Section 3 complete.\n")

# ============================================================================
# SECTION 4 — CROSS-COMPARISON
# ============================================================================

cat("\n=== SECTION 4: Cross-Comparison ===\n")

load_pio_sensitivity <- function(dir.sensitivity, ct_clean, model_label = "base") {
  fp <- paste0(dir.sensitivity, 'NEBULA/NEBULA_', ct_clean, '_', model_label, '.csv')
  if (!file.exists(fp)) { cat("  [load_pio] Not found:", fp, "\n"); return(NULL) }
  df <- read.csv(fp)
  df %>% transmute(gene   = summary.gene,
                   logFC  = summary.logFC_group_labelsYes,
                   pvalue = summary.p_group_labelsYes,
                   fdr    = p.adjust(pvalue, method = "BH"))
}

# ---- 4a. PIO vs GLP-1 scatter ----------------------------------------------
overlap_rows <- list()
for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  pio_df  <- load_pio_sensitivity(dir.sensitivity, ct_clean, "base")
  glp1_df <- load_glp1_effect(dir.root, ct_clean, "base")
  if (is.null(pio_df) || is.null(glp1_df)) next
  pio_df  <- pio_df  %>% filter(!is.na(logFC) & !is.na(pvalue))
  glp1_df <- glp1_df %>% filter(!is.na(logFC) & !is.na(pvalue))
  merged  <- inner_join(pio_df  %>% dplyr::select(gene, logFC_pio  = logFC, pval_pio  = pvalue),
                        glp1_df %>% dplyr::select(gene, logFC_glp1 = logFC, pval_glp1 = pvalue),
                        by = "gene") %>%
    mutate(sig_pio  = pval_pio  < P_THRESH,
           sig_glp1 = pval_glp1 < P_THRESH,
           status   = case_when(sig_pio & sig_glp1  ~ "Sig in both",
                                sig_pio & !sig_glp1 ~ "PIO only",
                                !sig_pio & sig_glp1 ~ "GLP-1 only",
                                TRUE                ~ "NS in both"),
           concordant = ifelse(sig_pio & sig_glp1, sign(logFC_pio) == sign(logFC_glp1), NA))
  if (nrow(merged) == 0) next
  r_val <- round(cor(merged$logFC_pio, merged$logFC_glp1, use = "complete.obs"), 3)
  overlap_rows[[ct]] <- data.frame(celltype = ct, n_genes_tested = nrow(merged),
                                   n_sig_pio = sum(merged$sig_pio), n_sig_glp1 = sum(merged$sig_glp1),
                                   n_sig_both = sum(merged$status == "Sig in both"),
                                   n_pio_only = sum(merged$status == "PIO only"),
                                   n_glp1_only = sum(merged$status == "GLP-1 only"),
                                   n_concordant = sum(merged$concordant, na.rm = TRUE),
                                   n_discordant = sum(!merged$concordant, na.rm = TRUE), r_logFC = r_val)
  delta   <- abs(merged$logFC_pio - merged$logFC_glp1)
  label_v <- ifelse(rank(-delta) <= 8 & merged$status != "NS in both", merged$gene, NA)
  p <- ggplot(merged, aes(x = logFC_pio, y = logFC_glp1, color = status, label = label_v)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggrepel::geom_text_repel(size = 2.5, color = "black", max.overlaps = 15, na.rm = TRUE) +
    scale_color_manual(values = c("Sig in both" = "#d73027", "PIO only" = "#E69F00",
                                  "GLP-1 only" = "#56B4E9", "NS in both" = "grey80")) +
    theme_bw() +
    labs(x = "PIO LogFC", y = "GLP-1 LogFC", color = NULL,
         title = paste0(ct, ": PIO vs GLP-1"), subtitle = paste0("r = ", r_val)) +
    theme(plot.title = element_text(face = "bold", size = 11), legend.position = "bottom")
  ggsave(paste0(dir.root, 'Comparison/Scatter_PIOvsGLP1_', ct_clean, '.pdf'), p, width = 8, height = 7)
  ggsave(paste0(dir.root, 'Comparison/Scatter_PIOvsGLP1_', ct_clean, '.png'), p, width = 8, height = 7, dpi = 200)
  cat("Saved scatter for", ct, "\n")
}
overlap_summary <- bind_rows(overlap_rows)
if (nrow(overlap_summary) > 0)
  write.csv(overlap_summary, paste0(dir.root, 'Comparison/PIO_vs_GLP1_DEG_Overlap.csv'), row.names = FALSE)

# ---- 4b. Shared concordant DEGs --------------------------------------------
shared_degs <- list()
for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  pio_df  <- load_pio_sensitivity(dir.sensitivity, ct_clean, "base")
  glp1_df <- load_glp1_effect(dir.root, ct_clean, "base")
  if (is.null(pio_df) || is.null(glp1_df)) next
  shared <- inner_join(
    pio_df  %>% filter(!is.na(pvalue) & pvalue < P_THRESH) %>%
      dplyr::select(gene, logFC_pio = logFC, pval_pio = pvalue, fdr_pio = fdr),
    glp1_df %>% filter(!is.na(pvalue) & pvalue < P_THRESH) %>%
      dplyr::select(gene, logFC_glp1 = logFC, pval_glp1 = pvalue, fdr_glp1 = fdr),
    by = "gene") %>%
    filter(sign(logFC_pio) == sign(logFC_glp1)) %>%
    mutate(direction = ifelse(logFC_pio > 0, "Up in both", "Down in both"), celltype = ct) %>%
    arrange(pval_pio)
  if (nrow(shared) > 0) shared_degs[[ct]] <- shared
}
if (length(shared_degs) > 0) {
  shared_combined <- bind_rows(shared_degs)
  write.csv(shared_combined, paste0(dir.root, 'Comparison/SharedConcordantDEGs_PIO_GLP1.csv'), row.names = FALSE)
  cat("Concordant DEGs:\n"); print(shared_combined %>% count(celltype, direction))
}

# ---- 4c. Sig gene count summary -------------------------------------------
sig_count_rows <- list()
for (ct in celltypes_vec) {
  ct_clean <- str_replace_all(str_replace_all(ct, "/", "_"), "-", "_")
  pio_df  <- load_pio_sensitivity(dir.sensitivity, ct_clean, "base")
  glp1_df <- load_glp1_effect(dir.root, ct_clean, "base")
  int_fp  <- paste0(dir.interaction, 'NEBULA/Interaction_NEBULA_', ct_clean, '.csv')
  int_df  <- if (file.exists(int_fp)) {
    df <- read.csv(int_fp)
    logfc_c <- grep("logFC.*group_labels.*glp1ra|logFC.*glp1ra.*group_labels",
                    names(df), value = TRUE, ignore.case = TRUE)[1]
    p_c     <- grep("summary\\.p_.*group_labels.*glp1ra|summary\\.p_.*glp1ra.*group_labels",
                    names(df), value = TRUE, ignore.case = TRUE)[1]
    gene_c  <- if ("summary.gene" %in% names(df)) "summary.gene" else "gene"
    if (!is.na(logfc_c) && !is.na(p_c))
      df %>% transmute(gene = .data[[gene_c]], logFC = .data[[logfc_c]],
                       pvalue = .data[[p_c]], fdr = p.adjust(pvalue, "BH"))
    else NULL
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
cat("\nSig gene count summary:\n"); print(sig_count_df)

cat("\n=== GLP-1 Extended Analysis Complete ===\n")
cat("Root output directory:", dir.root, "\n")