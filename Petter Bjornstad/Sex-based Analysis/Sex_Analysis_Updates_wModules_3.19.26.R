################################################################################
# Sex-Based Factor / Module / Severity Analysis — Full Pipeline
#
# Sections:
#   1. NMF module discovery (pseudo-bulk → cell scoring)
#   2. Morpha-style literature module scoring
#   3. NEBULA with HbA1c + eGFR + BMI covariates (continuous)
#   4. Severity-stratified sex comparisons (quartile + PS-matching + NEBULA)
#   5. Model comparison: standard vs severity-adjusted LogFC
################################################################################

library(Seurat)
library(NMF)
library(nebula)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(ggrepel)
library(rstatix)
library(data.table)
library(readxl)
library(pheatmap)
library(viridis)
library(MatchIt)        # install.packages("MatchIt")
library(cobalt)         # install.packages("cobalt") — balance diagnostics

set.seed(123)
# ── Output directories ────────────────────────────────────────────────────────
dir_base   <- "C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/"
dir_module <- file.path(dir_base, "Module_Analysis/")
dir_sev    <- file.path(dir_base, "Severity_Analysis/")
for (d in c(dir_module, dir_sev))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ── Severity variables used throughout ───────────────────────────────────────
SEVERITY_VARS <- c("hba1c", "eGFR_CKD_epi", "bmi")

################################################################################
# DATA LOADING
################################################################################

load("C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData")
so_full <- so_kpmp_sc
remove(so_kpmp_sc)

so_full <- subset(so_full, subset = record_id != "CRC-55")
so_full <- subset(so_full, subset = group %in% c("Type_2_Diabetes", "Lean_Control"))

so_full$celltype2 <- case_when(
  so_full$KPMP_celltype %in% c("aPT", "PT-S1/S2", "PT-S3") ~ "PT",
  grepl("TAL", so_full$KPMP_celltype)                        ~ "TAL",
  grepl("EC-",  so_full$KPMP_celltype)                       ~ "EC",
  TRUE                                                        ~ so_full$KPMP_celltype
)

# ── Clinical data ─────────────────────────────────────────────────────────────
harmonized_data <- read.csv(
  "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv",
  na = ""
)

clin <- harmonized_data %>%
  dplyr::select(-dob) %>%
  arrange(date_of_screen) %>%
  reframe(
    across(where(negate(is.numeric)), ~ last(na.omit(.x))),
    across(where(is.numeric),         ~ mean(na.omit(.x), na.rm = TRUE)),
    .by = c(record_id, visit)
  ) %>%
  filter(visit == "baseline") %>%
  dplyr::select(record_id, sex, group, age, bmi, hba1c,
                eGFR_CKD_epi, acr_u,
                any_of(c("gfr_raw_plasma", "gfr_bsa_plasma", "sbp", "dbp")))

cat("Clinical data loaded:", nrow(clin), "participants\n")
cat("Severity completeness:\n")
print(clin %>% summarise(across(all_of(SEVERITY_VARS), ~ sum(!is.na(.)))))

################################################################################
# HELPER: subset Seurat object by broad cell type label
################################################################################
subset_celltype <- function(so, label) {
  if (label == "All") return(so)
  if (label == "IC")  return(subset(so, KPMP_celltype %in%
                                      c("cDC","cycT","CD4+ T","CD8+ T","NK","B","MON","MAC","MC")))
  subset(so, subset = celltype2 == label)
}

################################################################################
# SECTION 1 — NMF MODULE DISCOVERY
################################################################################
# Run NMF on pseudo-bulk (participant × gene) so modules are stable and
# not dominated by cell-count imbalance between T2D and LC.
# Then score individual cells with AddModuleScore.

run_nmf_modules <- function(so_obj, celltype_label, n_factors = 8,
                            n_hvg = 3000, outdir = dir_module) {
  
  cat("\n===== NMF Module Discovery:", celltype_label, "=====\n")
  so_ct <- subset_celltype(so_obj, celltype_label)
  DefaultAssay(so_ct) <- "RNA"
  
  # ── 1a. Pseudo-bulk matrix ────────────────────────────────────────────────
  meta       <- so_ct@meta.data
  counts_raw <- GetAssayData(so_ct, assay = "RNA", layer = "counts")
  pids       <- unique(meta$record_id)
  
  pb_mat <- do.call(cbind, lapply(pids, function(pid) {
    idx <- which(meta$record_id == pid)
    Matrix::rowSums(counts_raw[, idx, drop = FALSE])
  }))
  colnames(pb_mat) <- pids
  
  # log-CPM normalise
  pb_log <- log1p(sweep(pb_mat, 2, colSums(pb_mat), "/") * 1e6)
  
  # Top HVGs across participants
  top_genes <- names(sort(apply(pb_log, 1, var), decreasing = TRUE))[
    1:min(n_hvg, nrow(pb_log))]
  pb_sub <- pb_log[top_genes, ]
  
  # ── 1b. NMF ──────────────────────────────────────────────────────────────
  set.seed(42)
  nmf_res      <- nmf(pb_sub, rank = n_factors, method = "brunet",
                      nrun = 10, .options = "v")
  W            <- basis(nmf_res)   # genes × k
  H            <- coef(nmf_res)    # k × participants
  factor_names <- paste0("NMF_F", seq_len(n_factors))
  colnames(W)  <- factor_names
  rownames(H)  <- factor_names
  
  # ── 1c. Top genes per factor (for AddModuleScore) ─────────────────────────
  top50 <- lapply(seq_len(n_factors), function(f)
    rownames(W)[order(W[, f], decreasing = TRUE)[1:50]])
  names(top50) <- factor_names
  
  map_dfr(factor_names, ~ data.frame(
    factor  = .x,
    gene    = top50[[.x]],
    loading = W[top50[[.x]], .x]
  )) %>%
    write.csv(file.path(outdir, paste0("NMF_top50genes_", celltype_label, ".csv")),
              row.names = FALSE)
  
  # ── 1d. Score individual cells ────────────────────────────────────────────
  so_ct <- AddModuleScore(so_ct, features = top50, name = "NMF_F", ctrl = 100)
  # Seurat appends 1-indexed integers: NMF_F1 … NMF_F8
  score_cols <- paste0("NMF_F", seq_len(n_factors))
  
  # ── 1e. Participant-level heatmap of factor weights ───────────────────────
  # Only annotate participants that exist in clin AND have non-NA sex/group
  ann_df <- clin %>%
    filter(record_id %in% pids, !is.na(sex), !is.na(group)) %>%
    dplyr::select(record_id, sex, group) %>%
    distinct(record_id, .keep_all = TRUE) %>%
    column_to_rownames("record_id")
  
  # Subset H to only those participants so dimensions match
  pids_ann <- intersect(pids, rownames(ann_df))
  H_sub    <- H[, pids_ann, drop = FALSE]
  
  if (ncol(H_sub) > 1) {
    png(file.path(outdir, paste0("NMF_heatmap_", celltype_label, ".png")),
        width = 1400, height = 900)
    pheatmap(t(H_sub),
             annotation_col = ann_df[pids_ann, , drop = FALSE],
             scale          = "row",
             color          = viridis(100),
             main           = paste("NMF Factor Weights —", celltype_label),
             fontsize       = 10)
    dev.off()
  } else {
    cat("Skipping heatmap for", celltype_label,
        "— too few annotated participants.\n")
  }
  
  # ── 1f. Sex tests per factor, within each group ───────────────────────────
  score_data <- so_ct@meta.data %>%
    dplyr::select(record_id, sex, group, all_of(score_cols)) %>%
    group_by(record_id, sex, group) %>%
    summarise(across(all_of(score_cols), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop")
  
  sex_tests <- map_dfr(score_cols, function(sc) {
    map_dfr(c("Type_2_Diabetes", "Lean_Control"), function(grp) {
      sub <- score_data %>% filter(group == grp, !is.na(sex))
      if (n_distinct(sub$sex) < 2 || nrow(sub) < 6) return(NULL)
      wt <- wilcox.test(as.formula(paste(sc, "~ sex")), data = sub)
      data.frame(factor = sc, group = grp, p_value = wt$p.value,
                 mean_male   = mean(sub[[sc]][sub$sex == "Male"],   na.rm = TRUE),
                 mean_female = mean(sub[[sc]][sub$sex == "Female"], na.rm = TRUE))
    })
  }) %>% mutate(p_fdr = p.adjust(p_value, "BH"))
  
  write.csv(sex_tests,
            file.path(outdir, paste0("NMF_sex_tests_", celltype_label, ".csv")),
            row.names = FALSE)
  
  # ── 1g. Violin plots for nominally significant factors ────────────────────
  sig_facs <- sex_tests %>% filter(p_value < 0.05) %>% pull(factor) %>% unique()
  if (length(sig_facs) > 0) {
    vplots <- lapply(sig_facs, function(fac) {
      ggplot(score_data %>% filter(!is.na(sex)),
             aes(x = sex, y = .data[[fac]], fill = sex)) +
        geom_violin(alpha = 0.5, trim = FALSE) +
        geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
        stat_compare_means(method = "wilcox.test", label = "p.format") +
        facet_wrap(~ group) +
        scale_fill_manual(values = c(Female = "#D55E00", Male = "#0072B2")) +
        labs(title = fac, x = NULL, y = "Module Score") +
        theme_classic(base_size = 11) + theme(legend.position = "none")
    })
    ggsave(file.path(outdir, paste0("NMF_sig_factors_", celltype_label, ".png")),
           wrap_plots(vplots, ncol = 2),
           width = 10, height = 4 * ceiling(length(vplots) / 2), dpi = 300)
  }
  
  cat("NMF complete for", celltype_label, "\n")
  return(list(nmf = nmf_res, W = W, H = H,
              modules = top50, scores = score_data, tests = sex_tests))
}

nmf_results <- list()
for (ct in c("PT", "TAL", "EC"))
  nmf_results[[ct]] <- run_nmf_modules(so_full, ct, n_factors = 8)

################################################################################
# SECTION 2 — MORPHA-STYLE LITERATURE MODULE SCORING
# Replace gene lists with exact sets from whichever Morpha/KPMP paper you read.
################################################################################

morpha_modules <- list(
  PT_healthy   = c("CUBN","LRP2","SLC5A2","SLC5A12","SLC13A3",
                   "SLC22A8","SLC22A6","SLC34A1","ALDOB","FBP1"),
  PT_injured   = c("VCAM1","HAVCR1","SPP1","CD44","LCN2",
                   "CXCL1","CXCL8","TIMP1","VIM","FN1"),
  PT_adaptive  = c("SLC7A5","SLC1A5","PCNA","MCM2","MKI67",
                   "TOP2A","HMGB2","TK1","TYMS"),
  aPT_sig      = c("VCAM1","SPP1","HAVCR1","LCN2","CD44",
                   "PROM1","CD24","ITGB6","SLC6A19"),
  mito_complex = c("NDUFB3","NDUFB8","NDUFS1","NDUFA4","COX4I1",
                   "COX7A2","ATP5F1B","UQCRC1","CYCS")
)

score_morpha_modules <- function(so_obj, modules, celltype_label,
                                 outdir = dir_module) {
  cat("\n===== Morpha Module Scoring:", celltype_label, "=====\n")
  so_ct <- subset_celltype(so_obj, celltype_label)
  
  # Filter to genes present in object
  mods_ok <- Filter(function(g) length(g) >= 3,
                    lapply(modules, function(g) intersect(g, rownames(so_ct))))
  cat(length(mods_ok), "modules with >= 3 genes present\n")
  
  so_ct <- AddModuleScore(so_ct, features = mods_ok, name = "morpha_", ctrl = 100)
  # Rename columns to module names for clarity
  raw_cols  <- paste0("morpha_", seq_along(mods_ok))
  nice_cols <- paste0("morpha_", names(mods_ok))
  so_ct@meta.data <- so_ct@meta.data %>%
    rename(setNames(raw_cols, nice_cols))
  
  # Participant-level scores (mean per participant per KPMP subtype)
  score_data <- so_ct@meta.data %>%
    dplyr::select(record_id, sex, group, KPMP_celltype, all_of(nice_cols)) %>%
    group_by(record_id, sex, group) %>%
    summarise(across(all_of(nice_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  write.csv(score_data,
            file.path(outdir, paste0("Morpha_scores_", celltype_label, ".csv")),
            row.names = FALSE)
  
  # Sex tests within each group
  sex_tests <- map_dfr(nice_cols, function(sc) {
    map_dfr(c("Type_2_Diabetes", "Lean_Control"), function(grp) {
      sub <- score_data %>% filter(group == grp, !is.na(sex))
      if (n_distinct(sub$sex) < 2 || nrow(sub) < 6) return(NULL)
      wt <- wilcox.test(as.formula(paste0("`", sc, "` ~ sex")), data = sub)
      data.frame(module = sc, group = grp, p_value = wt$p.value,
                 mean_male   = mean(sub[[sc]][sub$sex == "Male"],   na.rm = TRUE),
                 mean_female = mean(sub[[sc]][sub$sex == "Female"], na.rm = TRUE))
    })
  }) %>% mutate(p_fdr = p.adjust(p_value, "BH"))
  
  write.csv(sex_tests,
            file.path(outdir, paste0("Morpha_sex_tests_", celltype_label, ".csv")),
            row.names = FALSE)
  
  # Faceted violin grid: module × group
  pd_long <- score_data %>%
    pivot_longer(all_of(nice_cols), names_to = "module", values_to = "score") %>%
    mutate(module = str_remove(module, "morpha_"))
  
  p <- ggplot(pd_long, aes(x = sex, y = score, fill = sex)) +
    geom_violin(alpha = 0.5, trim = FALSE) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
    facet_grid(module ~ group, scales = "free_y") +
    scale_fill_manual(values = c(Female = "#D55E00", Male = "#0072B2")) +
    labs(title = paste("Morpha Module Scores —", celltype_label),
         x = NULL, y = "Module Score") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", strip.text = element_text(size = 8))
  
  ggsave(file.path(outdir, paste0("Morpha_scores_", celltype_label, ".png")),
         p, width = 10, height = 3 * length(nice_cols), dpi = 300)
  
  return(list(scores = score_data, tests = sex_tests))
}

morpha_results <- list()
for (ct in c("PT", "TAL"))
  morpha_results[[ct]] <- score_morpha_modules(so_full, morpha_modules, ct)

################################################################################
# SECTION 3 — NEBULA WITH SEVERITY COVARIATES (CONTINUOUS)
# Model: ~ sex + hba1c_z + egfr_z + bmi_z + sex:hba1c_z
# All severity variables are z-scored so coefficients are directly comparable.
################################################################################

run_nebula_severity <- function(so_obj, celltype_label, clin_df,
                                sev_vars    = SEVERITY_VARS,
                                group_filter = "Type_2_Diabetes",
                                outdir       = dir_sev) {
  
  cat("\n===== NEBULA + Severity:", celltype_label, "=====\n")
  so_ct <- subset_celltype(so_obj, celltype_label)
  so_ct <- subset(so_ct, subset = group == group_filter)
  DefaultAssay(so_ct) <- "RNA"
  
  counts_raw <- round(GetAssayData(so_ct, layer = "counts"))
  meta       <- so_ct@meta.data %>%
    left_join(clin_df %>% dplyr::select(record_id, all_of(sev_vars)),
              by = "record_id")
  
  # Drop cells missing sex or any severity variable
  keep <- complete.cases(meta[, c("sex", sev_vars)])
  cat("Cells with complete data:", sum(keep), "/", nrow(meta), "\n")
  meta       <- meta[keep, ]
  counts_raw <- counts_raw[, keep]
  
  if (n_distinct(meta$sex) < 2 || ncol(counts_raw) < 50) {
    cat("Insufficient data, skipping.\n"); return(NULL)
  }
  
  # Z-score severity variables
  sev_z <- paste0(sev_vars, "_z")
  for (i in seq_along(sev_vars))
    meta[[sev_z[i]]] <- as.numeric(scale(meta[[sev_vars[i]]]))
  
  # Model: sex + each severity_z + sex × hba1c (primary interaction of interest)
  formula_str <- paste(
    "~ sex +",
    paste(sev_z, collapse = " + "),
    "+ sex:hba1c_z"
  )
  pred_mat <- model.matrix(as.formula(formula_str), data = meta)
  
  data_g <- group_cell(count  = counts_raw, id  = meta$kit_id,
                       pred   = pred_mat,   offset = meta$pooled_offset)
  if (is.null(data_g))
    data_g <- list(count = counts_raw, id = meta$kit_id,
                   pred = pred_mat, library = meta$pooled_offset)
  
  result <- nebula(count = data_g$count, id = data_g$id,
                   pred  = data_g$pred,  offset = data_g$library,
                   ncore = 1, reml = TRUE, model = "NBLMM",
                   output_re = TRUE, covariance = TRUE)
  
  res_df <- as.data.frame(result) %>%
    mutate(num_cells  = nrow(meta),
           num_male   = sum(meta$sex == "Male"),
           num_female = sum(meta$sex == "Female"))
  
  ct_safe <- str_replace_all(celltype_label, "[/-]", "_")
  fout <- file.path(outdir,
                    paste0("NEBULA_severity_", ct_safe, "_", group_filter, ".csv"))
  write.csv(res_df, fout, row.names = FALSE)
  cat("Saved →", basename(fout), "\n")
  return(res_df)
}

severity_nebula <- list()
for (ct in c("All", "PT", "TAL", "EC", "POD"))
  severity_nebula[[ct]] <- run_nebula_severity(so_full, ct, clin)

################################################################################
# SECTION 4 — SEVERITY-STRATIFIED SEX COMPARISONS
# Three parallel strategies run in one pipeline:
#   4A — Quartile stratification of module scores
#   4B — Propensity-score matching (M vs F matched on all 3 severity vars)
#   4C — NEBULA on PS-matched participants only
################################################################################

# ── 4A: Quartile stratification of module/NMF scores ─────────────────────────

stratify_by_severity <- function(score_df, clin_df,
                                 score_cols,
                                 severity_var = "hba1c",
                                 group_filter = "Type_2_Diabetes",
                                 outdir       = dir_sev) {
  
  cat("\n--- Quartile stratification on", severity_var, "---\n")
  
  combined <- score_df %>%
    left_join(clin_df %>% dplyr::select(record_id, all_of(severity_var), sex, group),
              by = "record_id") %>%
    filter(group == group_filter, !is.na(sex), !is.na(.data[[severity_var]])) %>%
    mutate(sev_q = ntile(.data[[severity_var]], 4),
           sev_label = paste0(severity_var, " Q", sev_q))
  
  # Wilcoxon test per score per quartile
  qstats <- map_dfr(score_cols, function(sc) {
    combined %>%
      group_by(sev_label) %>%
      filter(n_distinct(sex) == 2) %>%
      summarise(
        score      = sc,
        n_male     = sum(sex == "Male"),
        n_female   = sum(sex == "Female"),
        mean_m     = mean(.data[[sc]][sex == "Male"],   na.rm = TRUE),
        mean_f     = mean(.data[[sc]][sex == "Female"], na.rm = TRUE),
        p_wilcox   = tryCatch(
          wilcox.test(.data[[sc]][sex == "Male"],
                      .data[[sc]][sex == "Female"])$p.value,
          error = function(e) NA_real_),
        .groups = "drop"
      )
  }) %>% mutate(p_fdr = p.adjust(p_wilcox, "BH"))
  
  write.csv(qstats,
            file.path(outdir, paste0("Quartile_", severity_var, "_sex_tests.csv")),
            row.names = FALSE)
  
  # Faceted violin grid
  pd <- combined %>%
    pivot_longer(all_of(score_cols), names_to = "module", values_to = "score")
  
  p <- ggplot(pd, aes(x = sex, y = score, fill = sex)) +
    geom_violin(alpha = 0.45, trim = FALSE) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 2.8) +
    facet_grid(module ~ sev_label, scales = "free_y") +
    scale_fill_manual(values = c(Female = "#D55E00", Male = "#0072B2")) +
    labs(title  = paste("Module Scores by Sex across", severity_var, "Quartiles"),
         x = NULL, y = "Score") +
    theme_bw(base_size = 9) +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 7),
          strip.text.y = element_text(size = 7))
  
  ggsave(file.path(outdir,
                   paste0("Quartile_", severity_var, "_module_sex.png")),
         p, width = 14, height = 3 * length(score_cols), dpi = 300)
  
  return(qstats)
}

# Run stratification for each severity variable and each cell type with scores
for (sev in SEVERITY_VARS) {
  for (ct in names(morpha_results)) {
    sc_cols <- names(morpha_results[[ct]]$scores) %>%
      keep(~ startsWith(.x, "morpha_"))
    stratify_by_severity(
      score_df     = morpha_results[[ct]]$scores,
      clin_df      = clin,
      score_cols   = sc_cols,
      severity_var = sev,
      group_filter = "Type_2_Diabetes"
    )
  }
}

# ── 4B: Propensity-score matching (all 3 severity vars simultaneously) ────────

match_on_severity <- function(clin_df,
                              match_vars   = SEVERITY_VARS,
                              group_filter = "Type_2_Diabetes",
                              ratio        = 1,
                              outdir       = dir_sev) {
  
  cat("\n===== PS Matching on:", paste(match_vars, collapse = ", "), "=====\n")
  
  df <- clin_df %>%
    filter(group == group_filter,
           !is.na(sex),
           complete.cases(across(all_of(match_vars)))) %>%
    mutate(sex_bin = as.integer(sex == "Male"))
  
  if (n_distinct(df$sex) < 2) { cat("Need both sexes.\n"); return(NULL) }
  
  formula_match <- as.formula(paste("sex_bin ~", paste(match_vars, collapse = " + ")))
  m_out <- matchit(formula_match, data = df,
                   method = "nearest", ratio = ratio, distance = "glm")
  
  matched_df <- match.data(m_out)
  cat("Matched N:", nrow(matched_df), "| Sex breakdown:\n")
  print(table(matched_df$sex))
  
  # Balance diagnostics
  bal <- bal.tab(m_out, thresholds = c(m = 0.1))
  print(bal)
  
  # Balance plot
  png(file.path(outdir, "PSMatch_balance_plot.png"), width = 800, height = 600)
  plot(bal)
  dev.off()
  
  # Before/after comparison for each severity variable
  bal_summary <- map_dfr(match_vars, function(sv) {
    all_m  <- mean(df[[sv]][df$sex == "Male"],   na.rm = TRUE)
    all_f  <- mean(df[[sv]][df$sex == "Female"], na.rm = TRUE)
    mat_m  <- mean(matched_df[[sv]][matched_df$sex == "Male"],   na.rm = TRUE)
    mat_f  <- mean(matched_df[[sv]][matched_df$sex == "Female"], na.rm = TRUE)
    data.frame(
      variable        = sv,
      mean_male_before   = all_m,  mean_female_before   = all_f,
      diff_before        = all_m - all_f,
      mean_male_after    = mat_m,  mean_female_after    = mat_f,
      diff_after         = mat_m  - mat_f
    )
  })
  print(bal_summary)
  write.csv(bal_summary, file.path(outdir, "PSMatch_balance_summary.csv"),
            row.names = FALSE)
  write.csv(matched_df,  file.path(outdir, "PSMatch_matched_participants.csv"),
            row.names = FALSE)
  
  # Visual check: severity distributions before/after
  sev_long_all <- df %>%
    pivot_longer(all_of(match_vars), names_to = "variable", values_to = "value") %>%
    mutate(sample = "Before matching")
  sev_long_mat <- matched_df %>%
    pivot_longer(all_of(match_vars), names_to = "variable", values_to = "value") %>%
    mutate(sample = "After matching")
  sev_long <- bind_rows(sev_long_all, sev_long_mat)
  
  p_bal <- ggplot(sev_long, aes(x = sex, y = value, fill = sex)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.4, size = 1.5) +
    facet_grid(variable ~ sample, scales = "free_y") +
    scale_fill_manual(values = c(Female = "#D55E00", Male = "#0072B2")) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
    labs(title = "Severity Balance Before and After PS Matching",
         x = NULL, y = "Value") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
  
  ggsave(file.path(outdir, "PSMatch_balance_boxplots.png"),
         p_bal, width = 8, height = 3 * length(match_vars), dpi = 300)
  
  return(list(matched_data = matched_df, match_obj = m_out,
              balance = bal_summary))
}

matched_res <- match_on_severity(clin, match_vars = SEVERITY_VARS)

# ── 4C: NEBULA on PS-matched participants ─────────────────────────────────────

run_nebula_matched <- function(so_obj, celltype_label, keep_ids,
                               outdir = dir_sev) {
  cat("\n===== NEBULA (PS-matched):", celltype_label, "=====\n")
  
  so_ct <- subset_celltype(so_obj, celltype_label)
  so_ct <- subset(so_ct, subset = record_id %in% keep_ids)
  DefaultAssay(so_ct) <- "RNA"
  
  counts_raw <- round(GetAssayData(so_ct, layer = "counts"))
  meta       <- so_ct@meta.data
  keep       <- complete.cases(meta[, "sex", drop = FALSE])
  meta       <- meta[keep, ]; counts_raw <- counts_raw[, keep]
  
  if (n_distinct(meta$sex) < 2 || ncol(counts_raw) < 50) {
    cat("Insufficient data, skipping.\n"); return(NULL)
  }
  
  pred_mat <- model.matrix(~ sex, data = meta)
  data_g   <- group_cell(count = counts_raw, id = meta$kit_id,
                         pred = pred_mat, offset = meta$pooled_offset)
  if (is.null(data_g))
    data_g <- list(count = counts_raw, id = meta$kit_id,
                   pred = pred_mat, library = meta$pooled_offset)
  
  result <- nebula(count = data_g$count, id = data_g$id,
                   pred  = data_g$pred,  offset = data_g$library,
                   ncore = 1, reml = TRUE, model = "NBLMM",
                   output_re = TRUE, covariance = TRUE)
  
  res_df  <- as.data.frame(result)
  ct_safe <- str_replace_all(celltype_label, "[/-]", "_")
  fout    <- file.path(outdir, paste0("NEBULA_PSmatched_", ct_safe, ".csv"))
  write.csv(res_df, fout, row.names = FALSE)
  cat("Saved →", basename(fout), "\n")
  return(res_df)
}

matched_nebula <- list()
if (!is.null(matched_res)) {
  matched_ids <- matched_res$matched_data$record_id
  for (ct in c("PT", "TAL", "EC"))
    matched_nebula[[ct]] <- run_nebula_matched(so_full, ct, matched_ids)
}

################################################################################
# SECTION 5 — THREE-WAY MODEL COMPARISON
# Compares gene-level logFC from:
#   (a) Standard model (sex only)
#   (b) Severity-adjusted model (sex + HbA1c + eGFR + BMI)
#   (c) PS-matched model (severity-balanced cohort)
################################################################################

compare_three_models <- function(std_file, sev_file, match_file,
                                 celltype_label, top_n = 30,
                                 outdir = dir_sev) {
  
  read_nebula <- function(f, suffix) {
    if (is.null(f) || !file.exists(f)) return(NULL)
    fread(f) %>%
      dplyr::select(gene = summary.gene,
                    !!paste0("logFC_", suffix) := summary.logFC_sexMale,
                    !!paste0("p_",    suffix) := summary.p_sexMale) %>%
      mutate(!!paste0("p_fdr_", suffix) := p.adjust(.data[[paste0("p_", suffix)]], "BH"))
  }
  
  std   <- read_nebula(std_file,   "std")
  sev   <- read_nebula(sev_file,   "sev")
  mat   <- read_nebula(match_file, "mat")
  
  # Join whatever models are available
  merged <- std
  if (!is.null(sev))  merged <- inner_join(merged, sev,  by = "gene")
  if (!is.null(mat))  merged <- inner_join(merged, mat,  by = "gene")
  merged <- drop_na(merged)
  
  # Correlations
  if (!is.null(sev))
    cat(celltype_label, "r(std, sev):",
        round(cor(merged$logFC_std, merged$logFC_sev, use = "complete.obs"), 3), "\n")
  if (!is.null(mat))
    cat(celltype_label, "r(std, mat):",
        round(cor(merged$logFC_std, merged$logFC_mat, use = "complete.obs"), 3), "\n")
  
  # Classify genes by which models call them significant
  merged <- merged %>%
    mutate(
      sig_std = !is.na(p_fdr_std) & p_fdr_std < 0.05,
      sig_sev = !is.null(sev)   & !is.na(p_fdr_sev) & p_fdr_sev < 0.05,
      sig_mat = !is.null(mat)   & !is.na(p_fdr_mat) & p_fdr_mat < 0.05,
      robust  = sig_std & sig_sev & sig_mat,
      status  = case_when(
        robust                        ~ "Robust (all 3 models)",
        sig_std & sig_sev & !sig_mat  ~ "Std + Sev only",
        sig_std & !sig_sev & sig_mat  ~ "Std + Matched only",
        sig_std & !sig_sev & !sig_mat ~ "Standard only",
        !sig_std & (sig_sev | sig_mat)~ "Severity/Matched only",
        TRUE                          ~ "Not significant"
      )
    )
  
  label_genes <- merged %>%
    filter(status != "Not significant") %>%
    slice_min(p_fdr_std, n = top_n, with_ties = FALSE) %>%
    pull(gene)
  
  status_colors <- c(
    "Robust (all 3 models)"     = "#0072B2",
    "Std + Sev only"            = "#56B4E9",
    "Std + Matched only"        = "#009E73",
    "Standard only"             = "#E69F00",
    "Severity/Matched only"     = "#CC79A7",
    "Not significant"           = "grey80"
  )
  
  # Scatter: standard vs severity-adjusted (if available)
  plot_list <- list()
  if (!is.null(sev)) {
    plot_list[["std_vs_sev"]] <- ggplot(merged,
                                        aes(x = logFC_std, y = logFC_sev, color = status)) +
      geom_point(alpha = 0.5, size = 1.5) +
      geom_text_repel(data = merged %>% filter(gene %in% label_genes),
                      aes(label = gene), size = 2.3, max.overlaps = 20) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
      scale_color_manual(values = status_colors) +
      labs(title = paste(celltype_label, "— Standard vs Severity-adjusted"),
           x = "LogFC (standard)", y = "LogFC (severity-adjusted)") +
      theme_classic(base_size = 11) +
      theme(legend.position = "bottom",
            legend.title = element_blank())
  }
  
  if (!is.null(mat)) {
    plot_list[["std_vs_mat"]] <- ggplot(merged,
                                        aes(x = logFC_std, y = logFC_mat, color = status)) +
      geom_point(alpha = 0.5, size = 1.5) +
      geom_text_repel(data = merged %>% filter(gene %in% label_genes),
                      aes(label = gene), size = 2.3, max.overlaps = 20) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
      scale_color_manual(values = status_colors) +
      labs(title = paste(celltype_label, "— Standard vs PS-matched"),
           x = "LogFC (standard)", y = "LogFC (PS-matched)") +
      theme_classic(base_size = 11) +
      theme(legend.position = "bottom",
            legend.title = element_blank())
  }
  
  if (length(plot_list) > 0) {
    ggsave(file.path(outdir, paste0("ModelComparison_", celltype_label, ".png")),
           wrap_plots(plot_list, ncol = length(plot_list)),
           width = 8 * length(plot_list), height = 7, dpi = 300)
  }
  
  write.csv(merged,
            file.path(outdir, paste0("ModelComparison_results_", celltype_label, ".csv")),
            row.names = FALSE)
  
  # Summary table
  cat("\nSignificance classification summary —", celltype_label, ":\n")
  print(table(merged$status))
  
  return(merged)
}

# Run comparison for each cell type
comparison_results <- list()
for (ct in c("PT", "TAL", "EC")) {
  ct_safe <- str_replace_all(ct, "[/-]", "_")
  comparison_results[[ct]] <- compare_three_models(
    std_file   = file.path(dir_base, "T2D_Only/",
                           paste0("Full_NEBULA_", ct, "_cells__T2D_pooledoffset.csv")),
    sev_file   = file.path(dir_sev,
                           paste0("NEBULA_severity_", ct_safe, "_Type_2_Diabetes.csv")),
    match_file = file.path(dir_sev,
                           paste0("NEBULA_PSmatched_", ct_safe, ".csv")),
    celltype_label = ct
  )
}

cat("\n===== Full pipeline complete =====\n")
cat("Outputs written to:\n  ", dir_module, "\n  ", dir_sev, "\n")






#Next steps
################################################################################
# Follow-up Analyses — Sex Differences
# Requires objects already created by the main pipeline:
#   nmf_results, morpha_results, so_full, clin, dir_module, dir_sev
#
# Sections:
#   A. NMF_F3 TAL characterisation (top genes + pathway enrichment)
#   B. aPT_sig TAL × Disease group interaction (LC vs T2D, by sex)
#   C. PT pseudotime vs Morpha module score correlations
#   D. NMF factor weights vs clinical severity variables
################################################################################

library(clusterProfiler)   # install.packages / BiocManager::install
library(org.Hs.eg.db)
library(enrichplot)
library(ggcorrplot)        # install.packages("ggcorrplot")
library(broom)

dir_followup <- file.path(dir_base, "Followup_Analyses/")
dir.create(dir_followup, showWarnings = FALSE, recursive = TRUE)

################################################################################
# A. CHARACTERISE NMF_F3 IN TAL
#    - ranked gene list from factor loadings (W matrix)
#    - ORA on top 100 genes (GO BP + KEGG)
#    - GSEA on full ranked list
#    - Dot plots saved to dir_followup
################################################################################

characterise_nmf_factor <- function(nmf_res_ct, factor_name = "NMF_F3",
                                    celltype_label = "TAL",
                                    top_n_ora = 100,
                                    outdir = dir_followup) {
  
  cat("\n===== NMF Factor Characterisation:", factor_name, "-", celltype_label, "=====\n")
  
  W <- nmf_res_ct$W                              # genes × factors
  if (!factor_name %in% colnames(W))
    stop(factor_name, " not found in W matrix columns: ", paste(colnames(W), collapse = ", "))
  
  # ── Ranked gene list (all genes, sorted by loading) ──────────────────────
  loadings   <- sort(W[, factor_name], decreasing = TRUE)
  gene_syms  <- names(loadings)
  
  # Convert gene symbols → Entrez IDs (required by clusterProfiler)
  id_map <- bitr(gene_syms, fromType = "SYMBOL",
                 toType   = c("ENTREZID", "SYMBOL"),
                 OrgDb    = org.Hs.eg.db)
  
  # Save top 50 genes with loadings
  top50_df <- data.frame(gene = gene_syms[1:50],
                         loading = loadings[1:50])
  write.csv(top50_df,
            file.path(outdir, paste0(factor_name, "_", celltype_label, "_top50genes.csv")),
            row.names = FALSE)
  cat("Top 10 genes in", factor_name, ":\n")
  print(head(top50_df, 10))
  
  # ── ORA: top N genes ─────────────────────────────────────────────────────
  top_syms    <- gene_syms[1:min(top_n_ora, length(gene_syms))]
  top_entrez  <- id_map$ENTREZID[id_map$SYMBOL %in% top_syms]
  bg_entrez   <- id_map$ENTREZID                                # background = all detected
  
  ora_go <- enrichGO(gene          = top_entrez,
                     universe      = bg_entrez,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2,
                     readable      = TRUE)
  
  ora_kegg <- enrichKEGG(gene          = top_entrez,
                         universe      = bg_entrez,
                         organism      = "hsa",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05)
  
  # Save tables
  write.csv(as.data.frame(ora_go),
            file.path(outdir, paste0(factor_name, "_", celltype_label, "_ORA_GO.csv")),
            row.names = FALSE)
  write.csv(as.data.frame(ora_kegg),
            file.path(outdir, paste0(factor_name, "_", celltype_label, "_ORA_KEGG.csv")),
            row.names = FALSE)
  
  # ── GSEA: full ranked list ────────────────────────────────────────────────
  ranked_vec <- setNames(loadings, gene_syms)
  ranked_entrez <- id_map %>%
    filter(SYMBOL %in% names(ranked_vec)) %>%
    mutate(loading = ranked_vec[SYMBOL]) %>%
    arrange(desc(loading)) %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    { setNames(.$loading, .$ENTREZID) }
  
  gsea_go <- gseGO(geneList      = ranked_entrez,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   minGSSize     = 10,
                   maxGSSize     = 500,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   verbose       = FALSE)
  
  write.csv(as.data.frame(gsea_go),
            file.path(outdir, paste0(factor_name, "_", celltype_label, "_GSEA_GO.csv")),
            row.names = FALSE)
  
  # ── Plots ─────────────────────────────────────────────────────────────────
  plots <- list()
  
  # Bar plot of top loadings
  plots$loadings <- ggplot(top50_df[1:20, ], aes(x = reorder(gene, loading), y = loading)) +
    geom_col(fill = "#0072B2", alpha = 0.8) +
    coord_flip() +
    labs(title = paste("Top 20 genes —", factor_name, "(", celltype_label, ")"),
         x = NULL, y = "NMF Loading") +
    theme_classic(base_size = 11)
  
  if (nrow(as.data.frame(ora_go)) > 0) {
    plots$ora_dot <- dotplot(ora_go, showCategory = 20, font.size = 9) +
      ggtitle(paste("GO BP ORA —", factor_name, "(", celltype_label, ")"))
  }
  
  if (nrow(as.data.frame(gsea_go)) > 0) {
    # Top 10 positive NES (high-loading gene enrichment)
    gsea_df <- as.data.frame(gsea_go) %>% arrange(desc(NES))
    plots$gsea_bar <- ggplot(head(gsea_df, 20),
                             aes(x = reorder(Description, NES), y = NES,
                                 fill = p.adjust)) +
      geom_col() +
      coord_flip() +
      scale_fill_gradient(low = "#D55E00", high = "grey80",
                          name = "FDR", limits = c(0, 0.05)) +
      labs(title = paste("GSEA GO BP —", factor_name, "(", celltype_label, ")"),
           x = NULL, y = "NES") +
      theme_classic(base_size = 10)
  }
  
  # Combined save
  if (length(plots) > 0) {
    ggsave(file.path(outdir, paste0(factor_name, "_", celltype_label, "_enrichment.png")),
           wrap_plots(plots, ncol = 1),
           width = 10, height = 6 * length(plots), dpi = 300)
  }
  
  cat("Done. ORA GO terms:", nrow(as.data.frame(ora_go)),
      "| GSEA GO terms:", nrow(as.data.frame(gsea_go)), "\n")
  
  return(list(top50 = top50_df, ora_go = ora_go, ora_kegg = ora_kegg, gsea_go = gsea_go))
}

# Run for NMF_F3 in TAL (primary interest from Image 5)
# Also run for other significant/interesting factors you want to explore
nmf_f3_tal <- characterise_nmf_factor(nmf_results[["TAL"]],
                                      factor_name    = "NMF_F3",
                                      celltype_label = "TAL")

# To run for any other factor, e.g. NMF_F1 in EC:
# nmf_f1_ec <- characterise_nmf_factor(nmf_results[["EC"]], "NMF_F1", "EC")

################################################################################
# B. aPT_sig TAL × DISEASE GROUP INTERACTION
#    Tests whether the sex difference in aPT_sig is:
#      (i)  present in LC only, T2D only, or both
#      (ii) significantly amplified in T2D vs LC (interaction term)
#    Output: faceted violin + stats table + interaction model
################################################################################

test_module_group_sex_interaction <- function(morpha_res, module_col = "morpha_aPT_sig",
                                              celltype_label = "TAL",
                                              clin_df = clin,
                                              outdir = dir_followup) {
  
  cat("\n===== aPT_sig × Group Interaction:", celltype_label, "=====\n")
  
  score_df <- morpha_res$scores %>%
    filter(!is.na(sex), group %in% c("Type_2_Diabetes", "Lean_Control")) %>%
    mutate(group_fac = factor(group,
                              levels = c("Lean_Control", "Type_2_Diabetes"),
                              labels = c("LC", "T2D")),
           sex_fac = factor(sex, levels = c("Female", "Male")))
  
  if (!module_col %in% names(score_df))
    stop(module_col, " not found. Available: ", paste(names(score_df), collapse = ", "))
  
  # ── (i) Within-group sex tests ────────────────────────────────────────────
  within_group <- score_df %>%
    group_by(group_fac) %>%
    summarise(
      n_male    = sum(sex == "Male"),
      n_female  = sum(sex == "Female"),
      mean_m    = mean(.data[[module_col]][sex == "Male"],   na.rm = TRUE),
      mean_f    = mean(.data[[module_col]][sex == "Female"], na.rm = TRUE),
      delta_mf  = mean_m - mean_f,
      p_wilcox  = wilcox.test(.data[[module_col]] ~ sex)$p.value,
      .groups   = "drop"
    ) %>%
    mutate(p_fdr = p.adjust(p_wilcox, "BH"))
  
  cat("Within-group sex differences:\n"); print(within_group)
  
  # ── (ii) Interaction: sex × group (linear model, participant-level) ───────
  # Also test whether the sex-difference *magnitude* is larger in T2D
  lm_int <- lm(as.formula(paste(module_col, "~ sex_fac * group_fac")),
               data = score_df)
  int_tidy <- tidy(lm_int) %>%
    mutate(p_fdr = p.adjust(p.value, "BH"))
  cat("\nInteraction model:\n"); print(int_tidy)
  
  write.csv(within_group,
            file.path(outdir, paste0("aPT_sig_within_group_", celltype_label, ".csv")),
            row.names = FALSE)
  write.csv(int_tidy,
            file.path(outdir, paste0("aPT_sig_interaction_model_", celltype_label, ".csv")),
            row.names = FALSE)
  
  # ── (iii) Violin: sex × group ─────────────────────────────────────────────
  p_violin <- ggplot(score_df,
                     aes(x = sex_fac, y = .data[[module_col]], fill = sex_fac)) +
    geom_violin(alpha = 0.5, trim = FALSE) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.08, alpha = 0.5, size = 1.5) +
    stat_compare_means(method = "wilcox.test", label = "p.format",
                       label.x = 1.3, size = 3.5) +
    facet_wrap(~ group_fac) +
    scale_fill_manual(values = c(Female = "#D55E00", Male = "#0072B2")) +
    labs(title  = paste("aPT_sig Score by Sex and Disease Group —", celltype_label),
         x = NULL, y = "aPT_sig Module Score") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  
  # ── (iv) Delta plot: visualise amplification ──────────────────────────────
  delta_df <- within_group %>%
    dplyr::select(group_fac, delta_mf, p_wilcox)
  
  p_delta <- ggplot(delta_df, aes(x = group_fac, y = delta_mf, fill = group_fac)) +
    geom_col(alpha = 0.75, width = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_text(aes(label = paste0("p=", round(p_wilcox, 3))),
              vjust = -0.5, size = 4) +
    scale_fill_manual(values = c(LC = "#999999", T2D = "#CC0000")) +
    labs(title  = paste("Sex Difference (Male − Female) in aPT_sig —", celltype_label),
         x = NULL, y = "Δ (Male − Female) Module Score") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")
  
  ggsave(file.path(outdir, paste0("aPT_sig_interaction_", celltype_label, ".png")),
         p_violin / p_delta + plot_layout(heights = c(3, 1)),
         width = 8, height = 9, dpi = 300)
  
  return(list(within_group = within_group, interaction_model = int_tidy))
}

apt_interaction_TAL <- test_module_group_sex_interaction(
  morpha_results[["TAL"]],
  module_col     = "morpha_aPT_sig",
  celltype_label = "TAL"
)

# Also run for PT (where the Morpha plots showed no significance — worth checking LC vs T2D)
apt_interaction_PT <- test_module_group_sex_interaction(
  morpha_results[["PT"]],
  module_col     = "morpha_aPT_sig",
  celltype_label = "PT"
)

################################################################################
# C. PT PSEUDOTIME vs MORPHA MODULE CORRELATIONS
#    Correlates per-cell pseudotime (from Slingshot stored in metadata) with
#    Morpha module scores, split by sex and group.
#    Adjust `pseudotime_col` to match the column name in your Seurat metadata.
################################################################################

correlate_pseudotime_modules <- function(so_obj, morpha_modules_list,
                                         celltype_label  = "PT",
                                         pseudotime_col  = "slingPseudotime_1",
                                         outdir          = dir_followup) {
  
  cat("\n===== Pseudotime × Morpha Modules:", celltype_label, "=====\n")
  
  so_ct <- subset_celltype(so_obj, celltype_label)
  
  # Check pseudotime column exists
  if (!pseudotime_col %in% colnames(so_ct@meta.data)) {
    cat("Available metadata columns containing 'sling' or 'pseudo':\n")
    print(grep("sling|pseudo|time", colnames(so_ct@meta.data),
               ignore.case = TRUE, value = TRUE))
    stop("Pseudotime column '", pseudotime_col, "' not found. ",
         "Update pseudotime_col argument to the correct column name above.")
  }
  
  # Score modules on cells
  mods_ok <- Filter(function(g) length(g) >= 3,
                    lapply(morpha_modules_list,
                           function(g) intersect(g, rownames(so_ct))))
  so_ct <- AddModuleScore(so_ct, features = mods_ok,
                          name = "pseudocor_", ctrl = 100)
  score_cols_raw  <- paste0("pseudocor_", seq_along(mods_ok))
  score_cols_nice <- paste0("m_", names(mods_ok))
  so_ct@meta.data <- so_ct@meta.data %>%
    rename(setNames(score_cols_raw, score_cols_nice))
  
  # Cell-level data frame
  cell_df <- so_ct@meta.data %>%
    dplyr::select(record_id, sex, group,
                  pseudotime = all_of(pseudotime_col),
                  all_of(score_cols_nice)) %>%
    filter(!is.na(pseudotime), !is.na(sex))
  
  cat("Cells with pseudotime:", nrow(cell_df), "\n")
  
  # ── Spearman correlations: pseudotime ~ each module, stratified by sex × group ──
  strata <- cell_df %>%
    distinct(sex, group) %>%
    arrange(sex, group)
  
  cor_results <- map_dfr(seq_len(nrow(strata)), function(i) {
    sx  <- strata$sex[i];   grp <- strata$group[i]
    sub <- cell_df %>% filter(sex == sx, group == grp)
    map_dfr(score_cols_nice, function(mc) {
      ct_res <- cor.test(sub$pseudotime, sub[[mc]],
                         method = "spearman", exact = FALSE)
      data.frame(sex = sx, group = grp, module = mc,
                 rho = ct_res$estimate, p_value = ct_res$p.value,
                 n   = nrow(sub))
    })
  }) %>% mutate(p_fdr = p.adjust(p_value, "BH"))
  
  write.csv(cor_results,
            file.path(outdir, paste0("Pseudotime_module_correlations_", celltype_label, ".csv")),
            row.names = FALSE)
  
  cat("Significant correlations (FDR < 0.05):\n")
  print(cor_results %>% filter(p_fdr < 0.05) %>% arrange(p_fdr))
  
  # ── Scatter plots for significant pairs ───────────────────────────────────
  sig_pairs <- cor_results %>% filter(p_fdr < 0.05)
  
  if (nrow(sig_pairs) > 0) {
    # One plot per module, faceted by sex × group
    sig_modules <- unique(sig_pairs$module)
    scatter_plots <- lapply(sig_modules, function(mc) {
      ggplot(cell_df, aes(x = pseudotime, y = .data[[mc]], color = sex)) +
        geom_point(alpha = 0.15, size = 0.6) +
        geom_smooth(method = "loess", se = TRUE, linewidth = 0.9) +
        facet_grid(sex ~ group) +
        scale_color_manual(values = c(Female = "#D55E00", Male = "#0072B2")) +
        labs(title = paste(mc, "vs Pseudotime —", celltype_label),
             x     = "Slingshot Pseudotime",
             y     = paste(mc, "Score")) +
        theme_bw(base_size = 10) +
        theme(legend.position = "none")
    })
    
    ggsave(file.path(outdir, paste0("Pseudotime_scatter_", celltype_label, ".png")),
           wrap_plots(scatter_plots, ncol = 2),
           width = 12, height = 5 * ceiling(length(scatter_plots) / 2), dpi = 300)
  } else {
    cat("No significant pseudotime–module correlations after FDR correction.\n")
  }
  
  # ── Correlation heatmap (all strata × modules) ────────────────────────────
  cor_wide <- cor_results %>%
    mutate(stratum = paste0(sex, "\n", group)) %>%
    dplyr::select(stratum, module, rho) %>%
    pivot_wider(names_from = module, values_from = rho) %>%
    column_to_rownames("stratum")
  
  p_wide <- cor_results %>%
    mutate(stratum = paste0(sex, "\n", group)) %>%
    dplyr::select(stratum, module, p_fdr) %>%
    pivot_wider(names_from = module, values_from = p_fdr) %>%
    column_to_rownames("stratum")
  
  png(file.path(outdir, paste0("Pseudotime_rho_heatmap_", celltype_label, ".png")),
      width = 900, height = 500)
  ggcorrplot(as.matrix(cor_wide),
             p.mat   = as.matrix(p_wide),
             sig.lvl = 0.05,
             insig   = "blank",
             lab     = TRUE,
             lab_size = 3,
             title   = paste("Pseudotime ~ Module Spearman ρ —", celltype_label)) %>%
    print()
  dev.off()
  
  return(list(correlations = cor_results, cell_data = cell_df))
}

# !! Update pseudotime_col to the actual column in your Seurat metadata !!
# Common names: "slingPseudotime_1", "Pseudotime", "PT_pseudotime"
pt_pseudotime_cor <- correlate_pseudotime_modules(
  so_full,
  morpha_modules_list = morpha_modules,
  celltype_label      = "PT",
  pseudotime_col      = "slingPseudotime_1"   # <-- adjust if needed
)

################################################################################
# D. NMF FACTOR WEIGHTS vs CLINICAL SEVERITY VARIABLES
#    For each cell type: correlate participant-level NMF scores with
#    HbA1c, eGFR, BMI — and test whether correlations differ by sex.
#    Output: correlation matrix heatmap + scatter plots for sig pairs
################################################################################

correlate_nmf_clinical <- function(nmf_res_ct, clin_df,
                                   celltype_label  = "TAL",
                                   sev_vars        = SEVERITY_VARS,
                                   group_filter    = "Type_2_Diabetes",
                                   outdir          = dir_followup) {
  
  cat("\n===== NMF × Clinical Variables:", celltype_label, "=====\n")
  
  # Participant-level NMF scores (from nmf_results)
  score_df <- nmf_res_ct$scores %>%
    filter(group == group_filter, !is.na(sex)) %>%
    left_join(clin_df %>%
                dplyr::select(record_id, all_of(sev_vars)),
              by = "record_id")
  
  factor_cols <- names(score_df) %>% keep(~ startsWith(.x, "NMF_F"))
  cat("Factors:", paste(factor_cols, collapse = ", "), "\n")
  cat("N participants:", nrow(score_df), "\n")
  
  # ── Spearman correlations: each NMF factor ~ each clinical variable ───────
  cor_full <- map_dfr(factor_cols, function(fc) {
    map_dfr(sev_vars, function(sv) {
      complete <- score_df %>% filter(!is.na(.data[[fc]]), !is.na(.data[[sv]]))
      if (nrow(complete) < 6) return(NULL)
      ct_res <- cor.test(complete[[fc]], complete[[sv]],
                         method = "spearman", exact = FALSE)
      data.frame(factor = fc, clinical_var = sv,
                 group  = "Combined",
                 sex    = "All",
                 rho    = ct_res$estimate, p_value = ct_res$p.value,
                 n      = nrow(complete))
    })
  })
  
  # Stratify by sex
  cor_sex <- map_dfr(c("Female", "Male"), function(sx) {
    sub <- score_df %>% filter(sex == sx)
    map_dfr(factor_cols, function(fc) {
      map_dfr(sev_vars, function(sv) {
        complete <- sub %>% filter(!is.na(.data[[fc]]), !is.na(.data[[sv]]))
        if (nrow(complete) < 5) return(NULL)
        ct_res <- cor.test(complete[[fc]], complete[[sv]],
                           method = "spearman", exact = FALSE)
        data.frame(factor = fc, clinical_var = sv,
                   group  = group_filter, sex = sx,
                   rho    = ct_res$estimate, p_value = ct_res$p.value,
                   n      = nrow(complete))
      })
    })
  })
  
  all_cor <- bind_rows(cor_full, cor_sex) %>%
    mutate(p_fdr = p.adjust(p_value, "BH"))
  
  write.csv(all_cor,
            file.path(outdir,
                      paste0("NMF_clinical_correlations_", celltype_label, ".csv")),
            row.names = FALSE)
  
  cat("Significant correlations (FDR < 0.05):\n")
  print(all_cor %>% filter(p_fdr < 0.05) %>% arrange(p_fdr))
  
  # ── Heatmap: rho matrix (all participants) ─────────────────────────────────
  rho_mat <- cor_full %>%
    dplyr::select(factor, clinical_var, rho) %>%
    pivot_wider(names_from = clinical_var, values_from = rho) %>%
    column_to_rownames("factor") %>%
    as.matrix()
  
  p_mat <- cor_full %>%
    left_join(all_cor %>% filter(sex == "All") %>%
                dplyr::select(factor, clinical_var, p_fdr),
              by = c("factor", "clinical_var")) %>%
    dplyr::select(factor, clinical_var, p_fdr) %>%
    pivot_wider(names_from = clinical_var, values_from = p_fdr) %>%
    column_to_rownames("factor") %>%
    as.matrix()
  
  p_heat <- ggcorrplot(rho_mat,
                       p.mat    = p_mat,
                       sig.lvl  = 0.05,
                       insig    = "blank",
                       lab      = TRUE, lab_size = 3.5,
                       colors   = c("#D55E00", "white", "#0072B2"),
                       title    = paste("NMF Factor × Clinical Severity ρ —",
                                        celltype_label))
  ggsave(file.path(outdir,
                   paste0("NMF_clinical_rho_heatmap_", celltype_label, ".png")),
         p_heat, width = max(5, length(sev_vars) * 2),
         height = max(4, length(factor_cols) * 0.6 + 2), dpi = 300)
  
  # ── Sex-stratified comparison for significant pairs ────────────────────────
  sig_pairs <- all_cor %>%
    filter(sex == "All", p_fdr < 0.05) %>%
    dplyr::select(factor, clinical_var)
  
  if (nrow(sig_pairs) > 0) {
    scatter_plots <- lapply(seq_len(nrow(sig_pairs)), function(i) {
      fc <- sig_pairs$factor[i];  sv <- sig_pairs$clinical_var[i]
      sub <- score_df %>% filter(!is.na(.data[[fc]]), !is.na(.data[[sv]]))
      ggplot(sub, aes(x = .data[[sv]], y = .data[[fc]], color = sex)) +
        geom_point(alpha = 0.7, size = 2.5) +
        geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
        scale_color_manual(values = c(Female = "#D55E00", Male = "#0072B2")) +
        stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label..,
                                                        sep = "~`,`~")),
                 size = 3) +
        labs(title = paste(fc, "~", sv, "|", celltype_label),
             x = sv, y = paste(fc, "Score")) +
        theme_classic(base_size = 11)
    })
    
    ggsave(file.path(outdir,
                     paste0("NMF_clinical_scatterplots_", celltype_label, ".png")),
           wrap_plots(scatter_plots, ncol = 3),
           width = 15, height = 5 * ceiling(length(scatter_plots) / 3), dpi = 300)
  }
  
  return(all_cor)
}

# Run for all three cell types
nmf_clinical_cor <- list()
for (ct in c("PT", "TAL", "EC"))
  nmf_clinical_cor[[ct]] <- correlate_nmf_clinical(nmf_results[[ct]],
                                                   clin_df        = clin,
                                                   celltype_label = ct)

# ── Combined summary across all cell types ────────────────────────────────────
summary_all <- bind_rows(
  lapply(names(nmf_clinical_cor), function(ct)
    nmf_clinical_cor[[ct]] %>% mutate(celltype = ct))
) %>%
  filter(p_fdr < 0.05) %>%
  arrange(p_fdr)

write.csv(summary_all,
          file.path(dir_followup, "NMF_clinical_correlations_SUMMARY.csv"),
          row.names = FALSE)

cat("\n===== Follow-up analyses complete =====\n")
cat("Outputs in:", dir_followup, "\n")
cat("\nKey questions to answer from outputs:\n")
cat("  A. What pathways define NMF_F3 in TAL? (check *_GSEA_GO.csv)\n")
cat("  B. Is aPT_sig sex difference amplified in T2D? (check *interaction_model*.csv)\n")
cat("  C. Does pseudotime correlate with injured/adaptive modules? (check *pseudotime_cor*.csv)\n")
cat("  D. Which NMF factors track with HbA1c/eGFR? (check *_SUMMARY.csv)\n")













# =============================================================================
# Sex Differences Deep-Dive Analysis
# Focus 1: aPT_sig Morpha module in TAL  (higher in males, p=0.038)
# Focus 2: NMF_F3 in TAL                (higher in females, p=0.0059)
#
# Builds on outputs from the existing sex_diff_analysis pipeline:
#   - NMF_top50genes_TAL.csv  (top NMF gene loadings)
#   - Morpha_scores_TAL.csv   (per-cell morpha scores)
# And directly from so_kpmp_sc for cell-level differential analysis.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(ggpubr)
  library(patchwork)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggrepel)
  library(viridis)
  library(broom)
  library(car)
  library(corrplot)
  library(RColorBrewer)
})

# =============================================================================
# PATHS — confirmed from previous sessions
# =============================================================================

seurat_rdata   <- "C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData"
dir_module     <- "C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/Module_Analysis/"
harmonized_csv <- paste0("C:/Users/netio/OneDrive - UW/Laura Pyle's files - ",
                         "Biostatistics Core Shared Drive/Data Harmonization/",
                         "Data Clean/harmonized_dataset.csv")
nebula_pt_csv  <- paste0("C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/",
                         "T2D_Only/Full_NEBULA_PT_cells__T2D_pooledoffset.csv")

out_dir <- file.path(dir_module, "sex_diff_deep_dive")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# HELPERS
# =============================================================================

save_plot <- function(p, fname, w = 8, h = 6) {
  ggsave(file.path(out_dir, fname), p, width = w, height = h, dpi = 150)
  message("  saved: ", fname)
}

# ORA helper — returns enrichResult or NULL
run_ora <- function(genes, bg_genes = NULL, ont = "BP") {
  ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = org.Hs.eg.db, drop = TRUE)$ENTREZID
  if (length(ids) < 5) { message("  Too few mappable genes for ORA"); return(NULL) }
  
  bg_ids <- NULL
  if (!is.null(bg_genes))
    bg_ids <- bitr(bg_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db, drop = TRUE)$ENTREZID
  
  enrichGO(gene          = ids,
           universe      = bg_ids,
           OrgDb         = org.Hs.eg.db,
           ont           = ont,
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.20,
           readable      = TRUE)
}

# KEGG ORA helper
run_kegg <- function(genes) {
  ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = org.Hs.eg.db, drop = TRUE)$ENTREZID
  if (length(ids) < 5) return(NULL)
  enrichKEGG(gene          = ids,
             organism      = "hsa",
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.20)
}

# =============================================================================
# 1. LOAD DATA
# =============================================================================

message("\n=== 1. Loading data ===")

# ── Seurat object
if (!exists("so_kpmp_sc")) {
  message("  Loading so_kpmp_sc from RData...")
  load(seurat_rdata)     # loads so_kpmp_sc
} else {
  message("  so_kpmp_sc already in environment")
}

# ── Clinical metadata (hba1c, bmi — egfr has 0 non-NA values in harmonized)
clin <- read_csv(harmonized_csv, show_col_types = FALSE) %>%
  reframe(
    across(where(negate(is.numeric)), ~ last(na.omit(.x))),
    across(where(is.numeric),         ~ mean(na.omit(.x), na.rm = TRUE)),
    .by = c(record_id, visit)
  ) %>%
  filter(visit == "baseline") %>%
  dplyr::select(record_id, sex, group, age, bmi, hba1c,
                eGFR_CKD_epi, acr_u,
                any_of(c("gfr_raw_plasma", "gfr_bsa_plasma", "sbp", "dbp"))) %>%
  distinct(record_id, .keep_all = TRUE)

# ── Morpha scores for TAL
morpha_tal_file <- file.path(dir_module, "Morpha_scores_TAL.csv")
stopifnot(file.exists(morpha_tal_file))
morpha_tal <- read_csv(morpha_tal_file, show_col_types = FALSE)
# Expected columns: cell_id, morpha_aPT_sig, morpha_PT_healthy, ..., sex, group, record_id
message("  Morpha TAL rows: ", nrow(morpha_tal), " | cols: ", paste(names(morpha_tal), collapse = ", "))

# ── NMF top genes for TAL
nmf_tal_file <- file.path(dir_module, "NMF_top50genes_TAL.csv")
stopifnot(file.exists(nmf_tal_file))
nmf_top_tal <- read_csv(nmf_tal_file, show_col_types = FALSE)
# Expected columns: factor, gene, loading
message("  NMF TAL top-gene file loaded. Factors: ", paste(unique(nmf_top_tal$factor), collapse = ", "))

# ── NMF factor scores (from Seurat metadata — aggregated per participant)
tal_cells <- so_kpmp_sc@meta.data %>%
  as_tibble(rownames = "cell_id") %>%
  filter(grepl("TAL", celltype2, ignore.case = TRUE)) %>%
  mutate(sex   = factor(sex,   levels = c("Female", "Male")),
         group = factor(group, levels = c("Lean_Control", "Type_2_Diabetes")))

nmf_cols_tal <- grep("^NMF_F", names(tal_cells), value = TRUE)
message("  TAL NMF columns found: ", paste(nmf_cols_tal, collapse = ", "))

# Sample-level NMF means for TAL
nmf_sample_tal <- tal_cells %>%
  group_by(record_id, sex, group) %>%
  summarise(across(all_of(nmf_cols_tal), ~ mean(.x, na.rm = TRUE),
                   .names = "{.col}"), .groups = "drop")

# =============================================================================
# 2. FOCUS 1: aPT_sig MODULE IN TAL — GENE-LEVEL DRIVERS
# =============================================================================

message("\n=== 2. aPT_sig TAL — gene-level driver analysis ===")

# ── 2A. Identify which genes belong to the aPT_sig module
# The morpha module gene list should be stored alongside the scores;
# try common naming conventions
apt_gene_file <- file.path(dir_module, "Morpha_aPT_sig_genes.csv")
apt_genes_module <- if (file.exists(apt_gene_file)) {
  read_csv(apt_gene_file, show_col_types = FALSE) %>% pull(1)
} else {
  # Fallback: use aPT_sig column genes from the Morpha object if saved,
  # or load from the Seurat assay if it exists as a gene set
  message("  aPT_sig gene list CSV not found — checking Seurat misc slot")
  tryCatch(
    so_kpmp_sc@misc$morpha_gene_sets$aPT_sig,
    error = function(e) {
      message("  Not in misc slot — please supply Morpha_aPT_sig_genes.csv")
      NULL
    }
  )
}

if (!is.null(apt_genes_module) && length(apt_genes_module) > 0) {
  message("  aPT_sig module gene count: ", length(apt_genes_module))
  
  # ── 2B. Pathway enrichment on the aPT_sig gene set (context: what is this module?)
  message("  Running GO + KEGG on aPT_sig gene set...")
  
  bg_genes <- rownames(so_kpmp_sc)  # use all detected genes as background
  
  go_apt_bp <- run_ora(apt_genes_module, bg_genes, ont = "BP")
  go_apt_cc <- run_ora(apt_genes_module, bg_genes, ont = "CC")
  kegg_apt  <- run_kegg(apt_genes_module)
  
  if (!is.null(go_apt_bp) && nrow(go_apt_bp) > 0) {
    p <- dotplot(go_apt_bp, showCategory = 20) +
      ggtitle("aPT_sig Gene Set — GO:BP Enrichment") +
      theme(axis.text.y = element_text(size = 8))
    save_plot(p, "aPT_sig_GO_BP.pdf", w = 10, h = 8)
    
    write_csv(as.data.frame(go_apt_bp),
              file.path(out_dir, "aPT_sig_GO_BP.csv"))
  }
  
  if (!is.null(kegg_apt) && nrow(kegg_apt) > 0) {
    p <- dotplot(kegg_apt, showCategory = 15) +
      ggtitle("aPT_sig Gene Set — KEGG Enrichment")
    save_plot(p, "aPT_sig_KEGG.pdf", w = 10, h = 7)
    write_csv(as.data.frame(kegg_apt),
              file.path(out_dir, "aPT_sig_KEGG.csv"))
  }
  
  # ── 2C. Within TAL cells: pseudo-bulk DE between M vs F for aPT_sig genes
  message("  Pseudo-bulk DE (Male vs Female) within TAL for aPT_sig genes...")
  
  # Subset TAL + genes in aPT_sig module
  tal_seurat <- subset(so_kpmp_sc, cells = tal_cells$cell_id)
  apt_genes_present <- intersect(apt_genes_module, rownames(tal_seurat))
  message("    aPT_sig genes present in TAL: ", length(apt_genes_present), "/",
          length(apt_genes_module))
  
  # Pseudo-bulk: sum counts per sample
  pb_mat <- AggregateExpression(
    tal_seurat,
    assays    = "RNA",
    slot      = "counts",
    group.by  = "record_id",
    return.seurat = FALSE
  )$RNA
  pb_mat <- pb_mat[apt_genes_present, , drop = FALSE]
  
  # Match sample metadata
  pb_meta <- tal_cells %>%
    distinct(record_id, sex, group) %>%
    filter(record_id %in% colnames(pb_mat)) %>%
    arrange(match(record_id, colnames(pb_mat)))
  pb_mat <- pb_mat[, pb_meta$record_id]
  
  # T2D only (where the sex difference was observed)
  t2d_idx  <- pb_meta$group == "Type_2_Diabetes"
  pb_t2d   <- pb_mat[, t2d_idx]
  meta_t2d <- pb_meta[t2d_idx, ]
  
  # limma-voom differential expression
  if (requireNamespace("limma", quietly = TRUE) &&
      requireNamespace("edgeR", quietly = TRUE)) {
    library(limma); library(edgeR)
    dge <- DGEList(counts = pb_t2d, group = meta_t2d$sex)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~ sex, data = meta_t2d)
    v   <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    de_apt <- topTable(fit, coef = "sexMale", n = Inf, sort.by = "P") %>%
      rownames_to_column("gene") %>%
      as_tibble()
    write_csv(de_apt, file.path(out_dir, "aPT_sig_TAL_DE_MvF_T2D.csv"))
    
    # Volcano: genes driving the male-higher aPT_sig score
    de_apt <- de_apt %>%
      mutate(sig = adj.P.Val < 0.05,
             label = ifelse(sig | abs(logFC) > 1, gene, ""))
    p_vol <- ggplot(de_apt, aes(logFC, -log10(P.Value),
                                color = sig, label = label)) +
      geom_point(alpha = 0.7, size = 1.8) +
      geom_text_repel(size = 2.5, max.overlaps = 20) +
      scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#D62728"),
                         labels = c("ns", "FDR<0.05")) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      labs(title = "aPT_sig genes in TAL: Male vs Female (T2D only)",
           subtitle = "Positive logFC = higher in males",
           x = "logFC (Male - Female)", y = "-log10(p)", color = NULL) +
      theme_bw(base_size = 11)
    save_plot(p_vol, "aPT_sig_TAL_volcano_MvF.pdf", w = 8, h = 6)
    
    # Heatmap of top DE genes
    top_apt_de <- de_apt %>% filter(adj.P.Val < 0.1) %>% arrange(P.Value)
    if (nrow(top_apt_de) >= 3) {
      mat_heat <- cpm(dge, log = TRUE)[top_apt_de$gene[1:min(40, nrow(top_apt_de))], ]
      ann_col  <- data.frame(Sex = meta_t2d$sex, row.names = colnames(mat_heat))
      pheatmap(mat_heat,
               annotation_col  = ann_col,
               scale           = "row",
               cluster_cols    = TRUE,
               show_colnames   = FALSE,
               color           = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
               main            = "Top DE genes in aPT_sig module\n(TAL cells, T2D only)",
               filename        = file.path(out_dir, "aPT_sig_TAL_heatmap_topDE.pdf"),
               width = 8, height = 7)
    }
    
    # Enrichment of male-higher genes within aPT_sig module
    male_up_apt <- de_apt %>% filter(logFC > 0, adj.P.Val < 0.2) %>% pull(gene)
    if (length(male_up_apt) >= 5) {
      go_male_up <- run_ora(male_up_apt, apt_genes_present, ont = "BP")
      if (!is.null(go_male_up) && nrow(go_male_up) > 0) {
        p <- dotplot(go_male_up, showCategory = 15) +
          ggtitle("Male-enriched aPT_sig genes in TAL — GO:BP")
        save_plot(p, "aPT_sig_TAL_male_enriched_GO.pdf", w = 10, h = 7)
        write_csv(as.data.frame(go_male_up),
                  file.path(out_dir, "aPT_sig_TAL_male_GO.csv"))
      }
    }
  } else {
    message("  limma/edgeR not available — skipping pseudo-bulk DE")
  }
  
} else {
  message("  aPT_sig gene list unavailable — skipping section 2B–C")
  message("  Supply 'Morpha_aPT_sig_genes.csv' in dir_module to enable")
}

# ── 2D. aPT_sig score: sex × T2D interaction + clinical correlation
message("  aPT_sig score: interaction model + clinical correlations...")

if ("morpha_aPT_sig" %in% names(morpha_tal)) {
  
  # Cell-level data — aggregate to sample means
  apt_sample <- morpha_tal %>%
    group_by(record_id, sex, group) %>%
    summarise(aPT_sig = mean(morpha_aPT_sig, na.rm = TRUE), .groups = "drop") %>%
    mutate(sex   = factor(sex,   levels = c("Female", "Male")),
           group = factor(group, levels = c("Lean_Control", "Type_2_Diabetes")))
  
  # Interaction model
  fit_int <- lm(aPT_sig ~ sex * group, data = apt_sample)
  sink(file.path(out_dir, "aPT_sig_TAL_lm_interaction.txt"))
  cat("=== aPT_sig TAL ~ sex * T2D ===\n")
  print(summary(fit_int))
  cat("\nType II ANOVA:\n")
  print(Anova(fit_int, type = "II"))
  sink()
  
  # 2×2 violin
  p_2x2 <- ggplot(apt_sample, aes(sex, aPT_sig, fill = sex)) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_jitter(width = 0.12, size = 1.5, alpha = 0.7) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    stat_compare_means(method = "wilcox.test", label = "p.format",
                       label.x = 1.4, size = 3.5) +
    facet_wrap(~ group, labeller = as_labeller(c(
      "Lean_Control"    = "Lean Control",
      "Type_2_Diabetes" = "Type 2 Diabetes"
    ))) +
    scale_fill_manual(values = c("Female" = "#E377C2", "Male" = "#1F77B4")) +
    labs(title   = "aPT_sig module score in TAL cells",
         subtitle = "Sex difference by disease group",
         x = NULL, y = "aPT_sig score") +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  save_plot(p_2x2, "aPT_sig_TAL_sex_x_group.pdf", w = 8, h = 5)
  
  # Clinical correlations for aPT_sig
  apt_clin <- apt_sample %>%
    left_join(clin, by = "record_id") %>%
    filter(group == "Type_2_Diabetes")
  
  clin_vars <- intersect(c("hba1c", "bmi", "age", "acr", "upcr",
                           "systolic_bp", "diastolic_bp"), names(apt_clin))
  message("    Clinical vars available for correlation: ", paste(clin_vars, collapse = ", "))
  
  corr_apt <- map_dfr(clin_vars, function(cv) {
    d <- apt_clin %>% filter(!is.na(.data[[cv]]))
    if (nrow(d) < 5) return(NULL)
    ct <- cor.test(d$aPT_sig, d[[cv]], method = "spearman", exact = FALSE)
    tibble(clinical_var = cv, rho = ct$estimate, p = ct$p.value,
           n = nrow(d))
  }) %>% mutate(fdr = p.adjust(p, method = "BH"))
  print(corr_apt)
  write_csv(corr_apt, file.path(out_dir, "aPT_sig_TAL_clinical_corr.csv"))
  
  # Scatter plots for significant correlations
  sig_clin_apt <- corr_apt %>% filter(p < 0.05) %>% pull(clinical_var)
  if (length(sig_clin_apt) > 0) {
    plots_apt <- map(sig_clin_apt, function(cv) {
      ggplot(apt_clin, aes(.data[[cv]], aPT_sig, color = sex)) +
        geom_point(alpha = 0.8, size = 2) +
        geom_smooth(method = "lm", se = TRUE, linetype = "dashed", linewidth = 0.8) +
        scale_color_manual(values = c("Female" = "#E377C2", "Male" = "#1F77B4")) +
        labs(title   = paste("aPT_sig TAL vs.", cv, "(T2D)"),
             subtitle = paste0("ρ=", round(corr_apt$rho[corr_apt$clinical_var == cv], 2),
                               ", p=", signif(corr_apt$p[corr_apt$clinical_var == cv], 2)),
             x = cv, y = "aPT_sig score") +
        theme_bw(base_size = 11)
    })
    p_clin_apt <- wrap_plots(plots_apt, ncol = 2)
    save_plot(p_clin_apt, "aPT_sig_TAL_clinical_scatter.pdf",
              w = 5 * min(2, length(sig_clin_apt)), h = 4 * ceiling(length(sig_clin_apt) / 2))
  }
  
  # ── NEBULA scRNA-seq cross-reference: do the top male-biased genes in aPT_sig
  #    also show male-biased expression in NEBULA PT analysis?
  if (file.exists(nebula_pt_csv) && exists("de_apt")) {
    message("  Cross-referencing with NEBULA PT sex differences...")
    nebula_pt <- read_csv(nebula_pt_csv, show_col_types = FALSE) %>%
      rename(gene_symbol  = summary.gene,
             logFC_nebula = summary.logFC_sexMale,
             p_nebula     = summary.p_sexMale) %>%
      mutate(fdr_nebula = p.adjust(p_nebula, method = "BH"))
    
    cross <- de_apt %>%
      rename(gene_symbol = gene) %>%
      left_join(nebula_pt %>% select(gene_symbol, logFC_nebula, p_nebula, fdr_nebula),
                by = "gene_symbol") %>%
      mutate(
        sig_TAL_DE = adj.P.Val < 0.2,
        sig_NEBULA = fdr_nebula < 0.05,
        consistent = sign(logFC) == sign(logFC_nebula),
        label = ifelse((sig_TAL_DE | sig_NEBULA) & consistent, gene_symbol, "")
      )
    
    p_cross <- ggplot(cross, aes(logFC_nebula, logFC,
                                 color = consistent & sig_TAL_DE,
                                 label = label)) +
      geom_point(alpha = 0.7, size = 1.8) +
      geom_text_repel(size = 2.2, max.overlaps = 25) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#D62728"),
                         labels = c("Inconsistent / ns", "Consistent & significant")) +
      labs(title    = "aPT_sig module genes: TAL pseudo-bulk vs. NEBULA PT",
           subtitle = "Positive = higher in males. Colored = consistent male-bias",
           x = "logFC in NEBULA PT (sexMale)", y = "logFC in TAL pseudo-bulk (Male−Female)",
           color = NULL) +
      theme_bw(base_size = 11)
    save_plot(p_cross, "aPT_sig_TAL_vs_NEBULA_PT.pdf", w = 8, h = 6)
    
    write_csv(cross, file.path(out_dir, "aPT_sig_cross_TAL_NEBULA.csv"))
    message("  Consistent male-biased genes: ",
            sum(cross$consistent & cross$sig_TAL_DE, na.rm = TRUE))
  }
}

# =============================================================================
# 3. FOCUS 2: NMF_F3 IN TAL — CHARACTERIZE THE FEMALE-ENRICHED PROGRAM
# =============================================================================

message("\n=== 3. NMF_F3 TAL — female-enriched program characterization ===")

f3_genes_raw <- nmf_top_tal %>%
  filter(grepl("F3|factor.*3|NMF_F3", factor, ignore.case = TRUE)) %>%
  arrange(desc(loading))
message("  F3 top genes rows: ", nrow(f3_genes_raw))
write_csv(f3_genes_raw, file.path(out_dir, "NMF_F3_TAL_top_genes.csv"))

f3_genes <- f3_genes_raw %>% pull(gene)
bg_genes  <- rownames(so_kpmp_sc)

# ── 3A. GO pathway enrichment (BP, MF, CC)
message("  Running GO enrichment for NMF_F3...")
go_f3_bp <- run_ora(f3_genes, bg_genes, ont = "BP")
go_f3_mf <- run_ora(f3_genes, bg_genes, ont = "MF")
kegg_f3  <- run_kegg(f3_genes)

for (res_obj in list(list(go_f3_bp, "NMF_F3_TAL_GO_BP", "NMF_F3 TAL — GO: Biological Process"),
                     list(go_f3_mf, "NMF_F3_TAL_GO_MF", "NMF_F3 TAL — GO: Molecular Function"),
                     list(kegg_f3,  "NMF_F3_TAL_KEGG",  "NMF_F3 TAL — KEGG Pathways"))) {
  obj   <- res_obj[[1]]
  fname <- res_obj[[2]]
  title <- res_obj[[3]]
  if (!is.null(obj) && nrow(obj) > 0) {
    p <- dotplot(obj, showCategory = 20) +
      ggtitle(title) +
      theme(axis.text.y = element_text(size = 8))
    save_plot(p, paste0(fname, ".pdf"), w = 11, h = 8)
    write_csv(as.data.frame(obj), file.path(out_dir, paste0(fname, ".csv")))
  }
}

# Cnetplot for top BP terms
if (!is.null(go_f3_bp) && nrow(go_f3_bp) > 0) {
  tryCatch({
    p_cnet <- cnetplot(go_f3_bp, showCategory = 8,
                       foldChange = setNames(f3_genes_raw$loading, f3_genes_raw$gene),
                       circular = FALSE, colorEdge = TRUE) +
      ggtitle("NMF_F3 TAL — Gene-Pathway Network (top 8 BP terms)")
    save_plot(p_cnet, "NMF_F3_TAL_cnetplot.pdf", w = 12, h = 10)
  }, error = function(e) message("  cnetplot skipped: ", e$message))
}

# ── 3B. Loading bar plot — top 30 genes in F3
p_load <- f3_genes_raw %>%
  slice_max(loading, n = 30) %>%
  mutate(gene = fct_reorder(gene, loading)) %>%
  ggplot(aes(loading, gene)) +
  geom_col(fill = "#2CA02C", alpha = 0.8) +
  labs(title    = "NMF_F3 TAL — Top 30 gene loadings",
       subtitle = "Female-enriched program (p=0.0059)",
       x = "NMF loading", y = NULL) +
  theme_bw(base_size = 11)
save_plot(p_load, "NMF_F3_TAL_top30_loadings.pdf", w = 7, h = 7)

# ── 3C. Heatmap of F3 top genes across samples
if ("NMF_F3" %in% nmf_cols_tal) {
  
  # Expression heatmap of top 30 F3 genes in TAL cells (pseudobulk)
  tal_seurat_obj <- subset(so_kpmp_sc, cells = tal_cells$cell_id)
  pb_f3 <- AggregateExpression(
    tal_seurat_obj,
    assays = "RNA", slot = "counts",
    group.by = "record_id",
    return.seurat = FALSE
  )$RNA
  
  f3_genes_present <- intersect(f3_genes[1:30], rownames(pb_f3))
  pb_f3_sub <- pb_f3[f3_genes_present, ]
  
  # Sample annotation
  ann_f3 <- tal_cells %>%
    distinct(record_id, sex, group) %>%
    filter(record_id %in% colnames(pb_f3_sub)) %>%
    left_join(nmf_sample_tal %>% select(record_id, NMF_F3), by = "record_id") %>%
    column_to_rownames("record_id") %>%
    mutate(NMF_F3_score = NMF_F3) %>%
    select(Sex = sex, Group = group, NMF_F3_score)
  
  ann_colors <- list(
    Sex   = c(Female = "#E377C2", Male = "#1F77B4"),
    Group = c(Lean_Control = "#2CA02C", Type_2_Diabetes = "#D62728")
  )
  
  pb_f3_ord <- pb_f3_sub[, rownames(ann_f3)]
  pheatmap(
    log1p(pb_f3_ord),
    annotation_col  = ann_f3,
    annotation_colors = ann_colors,
    scale           = "row",
    cluster_cols    = TRUE,
    show_colnames   = FALSE,
    color           = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    main            = "NMF_F3 Top Genes in TAL (pseudobulk)",
    filename        = file.path(out_dir, "NMF_F3_TAL_heatmap.pdf"),
    width = 10, height = 8
  )
  message("  F3 heatmap saved")
}

# ── 3D. Clinical correlations for NMF_F3 (sample-level means)
message("  NMF_F3 clinical correlations...")

if ("NMF_F3" %in% names(nmf_sample_tal)) {
  
  f3_clin <- nmf_sample_tal %>%
    select(record_id, sex, group, NMF_F3) %>%
    left_join(clin, by = "record_id")
  
  # All samples + T2D only
  for (grp_label in c("All", "T2D_only")) {
    d_sub <- if (grp_label == "T2D_only")
      filter(f3_clin, group == "Type_2_Diabetes") else f3_clin
    
    clin_vars <- intersect(c("hba1c", "bmi", "age", "acr", "upcr",
                             "systolic_bp", "diastolic_bp"), names(d_sub))
    
    corr_f3 <- map_dfr(clin_vars, function(cv) {
      d <- d_sub %>% filter(!is.na(.data[[cv]]))
      if (nrow(d) < 5) return(NULL)
      ct <- cor.test(d$NMF_F3, d[[cv]], method = "spearman", exact = FALSE)
      tibble(clinical_var = cv, rho = ct$estimate, p = ct$p.value, n = nrow(d))
    }) %>% mutate(fdr = p.adjust(p, method = "BH"),
                  subset = grp_label)
    print(corr_f3)
    write_csv(corr_f3,
              file.path(out_dir, paste0("NMF_F3_TAL_clinical_corr_", grp_label, ".csv")))
    
    # Bar plot of correlations
    if (nrow(corr_f3) > 0) {
      p_rho <- corr_f3 %>%
        mutate(clinical_var = fct_reorder(clinical_var, rho),
               sig = p < 0.05) %>%
        ggplot(aes(rho, clinical_var, fill = rho, alpha = sig)) +
        geom_col() +
        scale_fill_gradient2(low = "#1F77B4", high = "#D62728", mid = "white",
                             midpoint = 0, name = "ρ") +
        scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 1),
                           guide = "none") +
        geom_vline(xintercept = 0, linewidth = 0.5) +
        labs(title    = paste("NMF_F3 TAL — Clinical Correlations (Spearman)"),
             subtitle = paste0(grp_label, "  |  filled = p<0.05"),
             x = "ρ", y = NULL) +
        theme_bw(base_size = 11)
      save_plot(p_rho, paste0("NMF_F3_TAL_clinical_rho_", grp_label, ".pdf"), w = 6, h = 4)
    }
    
    # Scatter for significant terms
    sig_clin_f3 <- corr_f3 %>% filter(p < 0.05) %>% pull(clinical_var)
    if (length(sig_clin_f3) > 0) {
      plots_f3 <- map(sig_clin_f3, function(cv) {
        ggplot(d_sub, aes(.data[[cv]], NMF_F3, color = sex)) +
          geom_point(alpha = 0.8, size = 2) +
          geom_smooth(method = "lm", se = TRUE, linetype = "dashed", linewidth = 0.8) +
          scale_color_manual(values = c("Female" = "#E377C2", "Male" = "#1F77B4")) +
          labs(title   = paste("NMF_F3 TAL vs.", cv),
               subtitle = paste0("ρ=", round(corr_f3$rho[corr_f3$clinical_var == cv], 2),
                                 ", p=", signif(corr_f3$p[corr_f3$clinical_var == cv], 2)),
               x = cv, y = "NMF_F3 score") +
          theme_bw(base_size = 11)
      })
      p_f3_clin <- wrap_plots(plots_f3, ncol = 2)
      save_plot(p_f3_clin, paste0("NMF_F3_TAL_clinical_scatter_", grp_label, ".pdf"),
                w = 5 * min(2, length(sig_clin_f3)),
                h = 4 * ceiling(length(sig_clin_f3) / 2))
    }
  }
  
  # ── 3E. NMF_F3: sex × group violin
  p_f3_viol <- ggplot(f3_clin, aes(sex, NMF_F3, fill = sex)) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_jitter(width = 0.12, size = 1.5, alpha = 0.7) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    stat_compare_means(method = "wilcox.test", label = "p.format",
                       label.x = 1.4, size = 3.5) +
    facet_wrap(~ group, labeller = as_labeller(c(
      "Lean_Control"    = "Lean Control",
      "Type_2_Diabetes" = "Type 2 Diabetes"
    ))) +
    scale_fill_manual(values = c("Female" = "#E377C2", "Male" = "#1F77B4")) +
    labs(title    = "NMF_F3 score in TAL cells",
         subtitle = "Female-enriched program (overall p=0.0059)",
         x = NULL, y = "NMF_F3 score") +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  save_plot(p_f3_viol, "NMF_F3_TAL_sex_x_group.pdf", w = 8, h = 5)
  
  # ── 3F. Is NMF_F3 inversely correlated with aPT_sig? (they have opposite sex effects)
  if ("morpha_aPT_sig" %in% names(morpha_tal)) {
    apt_f3_compare <- morpha_tal %>%
      group_by(record_id, sex, group) %>%
      summarise(aPT_sig = mean(morpha_aPT_sig, na.rm = TRUE), .groups = "drop") %>%
      left_join(nmf_sample_tal %>% select(record_id, NMF_F3), by = "record_id")
    
    ct_compare <- cor.test(apt_f3_compare$aPT_sig, apt_f3_compare$NMF_F3,
                           method = "spearman", exact = FALSE)
    message("  aPT_sig vs NMF_F3 correlation: ρ=",
            round(ct_compare$estimate, 3), " p=", signif(ct_compare$p.value, 3))
    
    p_compare <- ggplot(apt_f3_compare, aes(NMF_F3, aPT_sig, color = sex, shape = group)) +
      geom_point(size = 2.5, alpha = 0.8) +
      geom_smooth(aes(group = 1), method = "lm", color = "grey40",
                  se = TRUE, linetype = "dashed") +
      scale_color_manual(values = c("Female" = "#E377C2", "Male" = "#1F77B4")) +
      scale_shape_manual(values = c("Lean_Control" = 16, "Type_2_Diabetes" = 17)) +
      annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
               label = paste0("ρ=", round(ct_compare$estimate, 2),
                              ", p=", signif(ct_compare$p.value, 3)),
               size = 3.5, color = "grey30") +
      labs(title    = "aPT_sig (male↑) vs. NMF_F3 (female↑) in TAL",
           subtitle = "Opposing sex-dimorphic programs",
           x = "NMF_F3 score (female-enriched)",
           y = "aPT_sig score (male-enriched)",
           color = "Sex", shape = "Group") +
      theme_bw(base_size = 11)
    save_plot(p_compare, "aPT_sig_vs_NMF_F3_TAL.pdf", w = 7, h = 5)
  }
}

# =============================================================================
# 4. SUMMARY TABLE
# =============================================================================

message("\n=== 4. Saving output summary ===")

summary_files <- list.files(out_dir, full.names = FALSE)
writeLines(
  c("Sex Differences Deep-Dive — Output Files",
    paste0("Generated: ", Sys.time()), "",
    "FOCUS 1: aPT_sig module in TAL (higher in males, p=0.038)",
    "  aPT_sig_GO_BP.csv/pdf             — Pathway enrichment of gene set",
    "  aPT_sig_KEGG.csv/pdf              — KEGG enrichment of gene set",
    "  aPT_sig_TAL_DE_MvF_T2D.csv        — Pseudo-bulk DE (Male vs Female, T2D only)",
    "  aPT_sig_TAL_volcano_MvF.pdf       — Volcano of within-module DE",
    "  aPT_sig_TAL_heatmap_topDE.pdf     — Heatmap of top DE genes",
    "  aPT_sig_TAL_male_enriched_GO.csv  — GO of male-higher aPT_sig genes",
    "  aPT_sig_TAL_lm_interaction.txt    — Sex × T2D interaction model",
    "  aPT_sig_TAL_sex_x_group.pdf       — 2×2 violin plot",
    "  aPT_sig_TAL_clinical_corr.csv     — Spearman vs clinical variables",
    "  aPT_sig_TAL_vs_NEBULA_PT.pdf      — Cross-reference with NEBULA PT sex DE",
    "",
    "FOCUS 2: NMF_F3 in TAL (higher in females, p=0.0059)",
    "  NMF_F3_TAL_top_genes.csv          — F3 gene loadings",
    "  NMF_F3_TAL_top30_loadings.pdf     — Bar chart of top loadings",
    "  NMF_F3_TAL_GO_BP/MF.csv/pdf       — GO enrichment",
    "  NMF_F3_TAL_KEGG.csv/pdf           — KEGG enrichment",
    "  NMF_F3_TAL_cnetplot.pdf           — Gene-pathway network",
    "  NMF_F3_TAL_heatmap.pdf            — Pseudobulk expression heatmap",
    "  NMF_F3_TAL_sex_x_group.pdf        — Violin by sex × group",
    "  NMF_F3_TAL_clinical_corr_*.csv    — Clinical correlations (All + T2D-only)",
    "  NMF_F3_TAL_clinical_rho_*.pdf     — Correlation bar charts",
    "  NMF_F3_TAL_clinical_scatter_*.pdf — Scatter plots",
    "  aPT_sig_vs_NMF_F3_TAL.pdf         — Inverse relationship plot"
  ),
  file.path(out_dir, "README.txt")
)

message("\n=== Done. All outputs in: ", out_dir, " ===")



