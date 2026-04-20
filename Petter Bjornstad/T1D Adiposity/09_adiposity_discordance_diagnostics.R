# T1D adiposity analysis -- BMI vs DXA discordance diagnostics
#
# Step 1 in the BMI/DXA discordance investigation.
# Goal: quantify how much of the BMI-vs-DXA difference in NEBULA results is
# driven by (a) ATTEMPT study composition and (b) measurement disagreement in
# subjects with both BMI and DXA.
#
# Produces:
#   results/diagnostics/study_by_obesity_subjects.csv
#   results/diagnostics/study_by_obesity_cells.csv
#   results/diagnostics/bmi_dxa_concordance_matched_subjects.csv
#   results/diagnostics/bmi_dxa_misclass_directionality.csv
#   results/diagnostics/figures/*
#
# Run after 02_meta_merge_w_seurat.R.

library(tidyverse)
library(dplyr)
library(ggplot2)
library(aws.s3)
library(jsonlite)
library(irr)        # kappa2()

user <- Sys.info()[["user"]]

if (user == "choiyej") {
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad/"
  keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
} else if (user == "yejichoi") {
  root_path <- "/mmfs1/gscratch/togo/yejichoi/"
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
} else if (user == "sleidholt") {
  root_path <- "/Users/Shared/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/"
  git_path <- "/Users/sleidholt/Documents/GitHub/CHCO-Code/Petter Bjornstad/T1D Adipostiy/"
} else {
  stop("Unknown user: please specify root path for this user.")
}

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

s3_bucket <- "t1d.adiposity"

# -----------------------------------------------------------------------------
# Load cell-level metadata (built by 03_scrna_visualizations.R) and subject-
# level clinical table.
# -----------------------------------------------------------------------------
meta_cells <- s3readRDS(
  object = "data_clean/pb90_attempt_subset_meta.rds",
  bucket = s3_bucket, region = ""
)

clin <- s3readRDS(
  object = "data_clean/t1d_hc_clinical_data.csv",
  bucket = s3_bucket, region = ""
)

# Subject-level table derived from the merged Seurat metadata.
subj <- meta_cells %>%
  distinct(record_id, .keep_all = TRUE) %>%
  dplyr::select(record_id, study, group, sex, age,
                bmi, bmip, dexa_body_fat,
                bmi_obesity, dxa_obesity)

cat(sprintf("Total subjects in scRNA set: %d\n", nrow(subj)))
cat("Study x group (subjects):\n"); print(table(subj$study, subj$group, useNA = "ifany"))

# -----------------------------------------------------------------------------
# A) Subject-level study x obesity cross-tabs
# -----------------------------------------------------------------------------
tab_subj_bmi <- subj %>%
  dplyr::count(study, group, bmi_obesity) %>%
  tidyr::pivot_wider(names_from = bmi_obesity, values_from = n, values_fill = 0) %>%
  dplyr::mutate(definition = "BMI") %>%
  dplyr::rename_with(~ paste0(.x, "_n"), any_of(c("Normal", "Overweight", "Obese")))

tab_subj_dxa <- subj %>%
  dplyr::count(study, group, dxa_obesity) %>%
  tidyr::pivot_wider(names_from = dxa_obesity, values_from = n, values_fill = 0) %>%
  dplyr::mutate(definition = "DXA") %>%
  dplyr::rename_with(~ paste0(.x, "_n"), any_of(c("Normal", "Overweight", "Obese")))

tab_subj <- bind_rows(tab_subj_bmi, tab_subj_dxa)
cat("\nSubject-level study x obesity:\n"); print(tab_subj, n = 40)

# -----------------------------------------------------------------------------
# B) Cell-level study x obesity cross-tabs (PT + overall)
# -----------------------------------------------------------------------------
cells_long <- meta_cells %>%
  dplyr::select(record_id, study, group, KPMP_celltype_general,
                bmi_obesity, dxa_obesity)

tab_cells_bmi <- cells_long %>%
  filter(!is.na(bmi_obesity)) %>%
  dplyr::count(KPMP_celltype_general, study, bmi_obesity, name = "n_cells") %>%
  dplyr::mutate(definition = "BMI") %>%
  dplyr::rename(obesity = bmi_obesity)

tab_cells_dxa <- cells_long %>%
  filter(!is.na(dxa_obesity)) %>%
  dplyr::count(KPMP_celltype_general, study, dxa_obesity, name = "n_cells") %>%
  dplyr::mutate(definition = "DXA") %>%
  dplyr::rename(obesity = dxa_obesity)

tab_cells <- bind_rows(tab_cells_bmi, tab_cells_dxa)

# % of cells in each obesity category contributed by each study
cell_share <- tab_cells %>%
  group_by(definition, KPMP_celltype_general, obesity) %>%
  dplyr::mutate(share = n_cells / sum(n_cells)) %>%
  ungroup()

cat("\nPT cell share by study within each obesity bucket (BMI & DXA):\n")
print(cell_share %>% filter(KPMP_celltype_general == "PT") %>%
        arrange(definition, obesity, desc(share)))

# -----------------------------------------------------------------------------
# C) Concordance / misclassification direction in the overlap subset
# -----------------------------------------------------------------------------
both <- subj %>% filter(!is.na(bmi_obesity) & !is.na(dxa_obesity))
cat(sprintf("\nSubjects with both BMI and DXA categories: %d\n", nrow(both)))

# Weighted kappa (ordinal)
k_unweighted <- irr::kappa2(both[, c("bmi_obesity", "dxa_obesity")])
k_weighted   <- irr::kappa2(both[, c("bmi_obesity", "dxa_obesity")], weight = "squared")
cat(sprintf("Cohen's kappa (unweighted): %.3f (p = %.3g)\n",
            k_unweighted$value, k_unweighted$p.value))
cat(sprintf("Cohen's kappa (quadratic-weighted): %.3f (p = %.3g)\n",
            k_weighted$value, k_weighted$p.value))

# 3x3 concordance matrix
concord_mat <- both %>%
  dplyr::count(bmi_obesity, dxa_obesity) %>%
  tidyr::pivot_wider(names_from = dxa_obesity, values_from = n,
                     values_fill = 0,
                     names_prefix = "dxa_")
cat("\nBMI (rows) x DXA (cols):\n"); print(concord_mat)

# Directional misclassification: is BMI under- or over-calling adiposity vs DXA?
misclass_dir <- both %>%
  dplyr::mutate(
    bmi_num = as.integer(bmi_obesity),
    dxa_num = as.integer(dxa_obesity),
    direction = case_when(
      bmi_num == dxa_num           ~ "Concordant",
      bmi_num <  dxa_num           ~ "BMI lower than DXA (BMI missing adiposity)",
      bmi_num >  dxa_num           ~ "BMI higher than DXA"
    )
  )

cat("\nMisclassification direction (overlap subjects):\n")
print(misclass_dir %>% dplyr::count(direction))

# Stratify direction by sex and study
cat("\nMisclassification direction by sex:\n")
print(misclass_dir %>% dplyr::count(sex, direction))
cat("\nMisclassification direction by study:\n")
print(misclass_dir %>% dplyr::count(study, direction))

# -----------------------------------------------------------------------------
# D) Cells per subject by study (to flag yield imbalance)
# -----------------------------------------------------------------------------
cells_per_subj <- meta_cells %>%
  dplyr::count(record_id, study, group, name = "n_cells")

yield_by_study <- cells_per_subj %>%
  group_by(study) %>%
  summarise(n_subjects = dplyr::n(),
            total_cells = sum(n_cells),
            median_cells_per_subj = median(n_cells),
            iqr_cells_per_subj = IQR(n_cells),
            .groups = "drop")
cat("\nYield by study:\n"); print(yield_by_study)

# -----------------------------------------------------------------------------
# E) Save diagnostics to S3
# -----------------------------------------------------------------------------
save_csv_s3 <- function(df, key) {
  tmp <- tempfile(fileext = ".csv")
  write.csv(df, tmp, row.names = FALSE)
  put_object(file = tmp, object = key, bucket = s3_bucket, region = "")
  unlink(tmp)
  cat(sprintf("Saved %s\n", key))
}

save_csv_s3(tab_subj,     "results/diagnostics/study_by_obesity_subjects.csv")
save_csv_s3(tab_cells,    "results/diagnostics/study_by_obesity_cells.csv")
save_csv_s3(cell_share,   "results/diagnostics/cell_share_by_study_per_obesity.csv")
save_csv_s3(concord_mat,  "results/diagnostics/bmi_dxa_concordance_matched_subjects.csv")
save_csv_s3(misclass_dir %>% dplyr::count(direction),
            "results/diagnostics/bmi_dxa_misclass_directionality.csv")
save_csv_s3(yield_by_study, "results/diagnostics/yield_by_study.csv")

kappa_summary <- data.frame(
  metric = c("n_overlap", "kappa_unweighted", "kappa_weighted_sq",
             "kappa_unweighted_pvalue", "kappa_weighted_pvalue"),
  value  = c(nrow(both), k_unweighted$value, k_weighted$value,
             k_unweighted$p.value, k_weighted$p.value)
)
save_csv_s3(kappa_summary, "results/diagnostics/bmi_dxa_kappa_summary.csv")

# -----------------------------------------------------------------------------
# F) Figures
# -----------------------------------------------------------------------------
upload_plot <- function(p, key, width = 9, height = 6) {
  tmp <- tempfile(fileext = ".png")
  ggsave(tmp, plot = p, width = width, height = height, dpi = 300)
  put_object(file = tmp, object = key, bucket = s3_bucket, region = "")
  unlink(tmp)
  cat(sprintf("Saved figure %s\n", key))
}

# (1) Stacked bar: cell share by study within each obesity bucket (PT)
p_pt_share <- cell_share %>%
  filter(KPMP_celltype_general == "PT") %>%
  ggplot(aes(x = obesity, y = share, fill = study)) +
  geom_col(position = "stack") +
  facet_wrap(~ definition) +
  labs(title = "PT cell share by study within each obesity category",
       y = "Fraction of cells", x = NULL, fill = "Study") +
  theme_minimal()
upload_plot(p_pt_share, "results/diagnostics/figures/pt_cell_share_by_study.png")

# (2) Heatmap of BMI x DXA concordance in overlap subjects
p_concord <- misclass_dir %>%
  dplyr::count(bmi_obesity, dxa_obesity) %>%
  ggplot(aes(x = dxa_obesity, y = bmi_obesity, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), size = 5) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
  labs(title = sprintf("BMI vs DXA classification (n=%d with both)", nrow(both)),
       subtitle = sprintf("Cohen's kappa = %.3f", k_unweighted$value),
       x = "DXA obesity", y = "BMI obesity") +
  theme_minimal()
upload_plot(p_concord, "results/diagnostics/figures/bmi_dxa_concordance_heatmap.png",
            width = 7, height = 6)

# (3) Cells per subject by study (log scale)
p_yield <- cells_per_subj %>%
  ggplot(aes(x = study, y = n_cells, color = study)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  scale_y_log10() +
  labs(title = "scRNA yield per subject by study",
       y = "Cells per subject (log10)", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none")
upload_plot(p_yield, "results/diagnostics/figures/yield_per_subject_by_study.png",
            width = 7, height = 5)

cat("\nDone. Review outputs in s3://t1d.adiposity/results/diagnostics/.\n")
