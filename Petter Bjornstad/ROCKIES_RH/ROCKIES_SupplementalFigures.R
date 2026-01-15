########################################################################
# ROCKIES Supplemental Figures
########################################################################

library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggbreak)
library(patchwork)
library(cowplot)
library(corrplot)
library(data.table)
library(gtsummary)
library(gt)
library(stringr)

base_path <- 'C:/Users/netio/Documents/UofW/Rockies/publication_figures/'

########################################################################
# DATA PREPARATION (include all 3 groups)
########################################################################

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))

PET_avg <- function(data){
  tmp_df <- data %>% dplyr::select(lc_k2, rc_k2, lm_k2, rm_k2, lc_f, rc_f, lm_f, rm_f)
  avg_c_k2 <- tmp_df %>% dplyr::select(lc_k2, rc_k2) %>% rowMeans(na.rm=T)
  avg_m_k2 <- tmp_df %>% dplyr::select(lm_k2, rm_k2) %>% rowMeans(na.rm=T)
  avg_c_f <- tmp_df %>% dplyr::select(lc_f, rc_f) %>% rowMeans(na.rm=T)
  avg_m_f <- tmp_df %>% dplyr::select(lm_f, rm_f) %>% rowMeans(na.rm=T)
  avg_c_k2_f <- avg_c_k2 / avg_c_f
  avg_m_k2_f <- avg_m_k2 / avg_m_f
  results <- bind_cols(avg_c_k2, avg_m_k2, avg_c_f, avg_m_f, avg_c_k2_f, avg_m_k2_f) %>% as.data.frame()
  names(results) <- c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 'avg_c_k2_f', 'avg_m_k2_f')
  return(results)
}

# Remove existing PET columns and recalculate
dat <- dat %>% dplyr::select(-any_of(c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 'avg_c_k2_f', 'avg_m_k2_f')))
tmp_results <- PET_avg(dat)
dat_all <- dat %>% bind_cols(tmp_results)

# Filter to those with PET data and all 3 groups
dat_all <- dat_all %>% filter(!is.na(avg_c_k2))
dat_all <- dat_all %>% filter(group %in% c('Lean Control', 'Obese Control', 'Type 2 Diabetes'))

# Fix any misclassified participants
dat_all$group[which(dat_all$record_id == 'RH2-39-O')] <- 'Obese Control'

# Remove specific outlier if needed
dat_all <- dat_all %>% filter(record_id != 'CRC-55')

########################################################################
# SUPPLEMENTAL FIGURE 1: PET Comparisons with All 3 Groups
########################################################################

# Prepare data for plotting
supp1_df <- dat_all %>% 
  dplyr::select(record_id, group, avg_c_f, avg_c_k2, avg_c_k2_f)

supp1_long <- supp1_df %>%
  pivot_longer(cols = c(avg_c_f, avg_c_k2, avg_c_k2_f), 
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, 
                         levels = c("avg_c_f", "avg_c_k2", "avg_c_k2_f"),
                         labels = c("Cortical F", "Cortical K2", "Cortical K2/F")),
         group = factor(group, levels = c("Lean Control", "Obese Control", "Type 2 Diabetes")))

# Calculate pairwise statistics
stat_test_supp1 <- supp1_long %>%
  group_by(metric) %>%
  pairwise_wilcox_test(value ~ group, p.adjust.method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "metric", dodge = 0.8)

# Create the plot
p_supp1 <- ggplot(supp1_long, aes(x = metric, y = value)) +
  geom_boxplot(aes(fill = group), color = "black", alpha = 0.8, outlier.shape = NA,
               linewidth = 0.5, position = position_dodge(width = 0.8)) +
  geom_point(aes(fill = group), 
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.15),
             alpha = 0.6, size = 2.5, shape = 21, color = "black", stroke = 0.3) +
  scale_fill_manual(values = c("Lean Control" = "#c2dfe3", 
                               "Obese Control" = "#e8d5b7", 
                               "Type 2 Diabetes" = "#fcb1a6"),
                    labels = c("Lean Control", "Obese Control", "Type 2 Diabetes")) +
  labs(x = "PET Metrics", y = "Value", fill = "Group",
       title = "Supplemental Figure 1.", 
       subtitle = "Global PET Imaging Metrics by Group (Including Obese Controls)") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 20, face = "bold", color = "black"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0, color = "black"),
    plot.subtitle = element_text(size = 14, hjust = 0, color = "black", margin = margin(b = 10)),
    legend.position = "bottom",
    legend.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16, face = "bold", color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Add statistics - adjust y positions for visibility
stat_test_supp1_adj <- stat_test_supp1 %>%
  mutate(y.position = case_when(
    metric == "Cortical K2" ~ y.position * 1.05,
    metric == "Cortical K2/F" ~ y.position * 1.05,
    TRUE ~ y.position
  ))

p_supp1_stats <- p_supp1 + 
  stat_pvalue_manual(stat_test_supp1_adj, 
                     label = "p.adj.signif", 
                     tip.length = 0.01, 
                     hide.ns = TRUE, 
                     size = 5,
                     step.increase = 0.08)

# Add break in y-axis
p_supp1_final <- p_supp1_stats + scale_y_break(c(0.35, 0.65), scales = 2)

ggsave(paste0(base_path, 'SuppFig1_PET_3Groups.pdf'), plot = p_supp1_final,
       width = 14, height = 12, device = cairo_pdf, dpi = 300)
ggsave(paste0(base_path, 'SuppFig1_PET_3Groups.png'), plot = p_supp1_final,
       width = 14, height = 12, dpi = 300)
ggsave(paste0(base_path, 'SuppFig1_PET_3Groups.tiff'), plot = p_supp1_final,
       width = 14, height = 12, dpi = 300, compression = "lzw")

print("Supplemental Figure 1 saved!")

# Print summary statistics
cat("\n=== SUPPLEMENTAL FIGURE 1: Summary Statistics ===\n\n")
summary_stats_supp1 <- supp1_long %>%
  group_by(metric, group) %>%
  summarise(
    n = n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    .groups = 'drop'
  )
print(summary_stats_supp1)

write.csv(summary_stats_supp1, 
          paste0(base_path, 'SuppFig1_summary_statistics.csv'), 
          row.names = FALSE)

########################################################################
# SUPPLEMENTAL FIGURE 2: Medullary PET Metrics
########################################################################

supp2_df <- dat_all %>% 
  dplyr::select(record_id, group, avg_m_f, avg_m_k2, avg_m_k2_f)

supp2_long <- supp2_df %>%
  pivot_longer(cols = c(avg_m_f, avg_m_k2, avg_m_k2_f), 
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, 
                         levels = c("avg_m_f", "avg_m_k2", "avg_m_k2_f"),
                         labels = c("Medullary F", "Medullary K2", "Medullary K2/F")),
         group = factor(group, levels = c("Lean Control", "Obese Control", "Type 2 Diabetes")))

# Remove NaN values
supp2_long <- supp2_long %>% filter(!is.nan(value) & !is.na(value))

# Calculate pairwise statistics
stat_test_supp2 <- supp2_long %>%
  group_by(metric) %>%
  pairwise_wilcox_test(value ~ group, p.adjust.method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "metric", dodge = 0.8)

p_supp2 <- ggplot(supp2_long, aes(x = metric, y = value)) +
  geom_boxplot(aes(fill = group), color = "black", alpha = 0.8, outlier.shape = NA,
               linewidth = 0.5, position = position_dodge(width = 0.8)) +
  geom_point(aes(fill = group), 
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.15),
             alpha = 0.6, size = 2.5, shape = 21, color = "black", stroke = 0.3) +
  scale_fill_manual(values = c("Lean Control" = "#c2dfe3", 
                               "Obese Control" = "#e8d5b7", 
                               "Type 2 Diabetes" = "#fcb1a6")) +
  labs(x = "PET Metrics", y = "Value", fill = "Group",
       title = "Supplemental Figure 2.", 
       subtitle = "Medullary PET Imaging Metrics by Group") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 20, face = "bold", color = "black"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0, color = "black"),
    plot.subtitle = element_text(size = 14, hjust = 0, color = "black", margin = margin(b = 10)),
    legend.position = "bottom",
    legend.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16, face = "bold", color = "black"),
    panel.background = element_rect(fill = "white", color = NA)
  )

p_supp2_stats <- p_supp2 + 
  stat_pvalue_manual(stat_test_supp2, 
                     label = "p.adj.signif", 
                     tip.length = 0.01, 
                     hide.ns = TRUE, 
                     size = 5,
                     step.increase = 0.08)

ggsave(paste0(base_path, 'SuppFig2_Medullary_PET.pdf'), plot = p_supp2_stats,
       width = 14, height = 10, device = cairo_pdf, dpi = 300)
ggsave(paste0(base_path, 'SuppFig2_Medullary_PET.png'), plot = p_supp2_stats,
       width = 14, height = 10, dpi = 300)
ggsave(paste0(base_path, 'SuppFig2_Medullary_PET.tiff'), plot = p_supp2_stats,
       width = 14, height = 10, dpi = 300, compression = "lzw")

print("Supplemental Figure 2 saved!")

########################################################################
# SUPPLEMENTAL FIGURE 3: PET Metrics Correlation Matrix
########################################################################

# Prepare correlation data
corr_df <- dat_all %>%
  dplyr::select(avg_c_k2, avg_c_f, avg_c_k2_f, avg_m_k2, avg_m_f, avg_m_k2_f)

# Rename for cleaner labels
names(corr_df) <- c("Cortical K2", "Cortical F", "Cortical K2/F", 
                    "Medullary K2", "Medullary F", "Medullary K2/F")

# Remove rows with NaN
corr_df <- corr_df %>% filter(complete.cases(.))

# Calculate correlation matrix
corr_matrix <- cor(corr_df, method = "spearman", use = "pairwise.complete.obs")

# Calculate p-values
cor.mtest <- function(mat, method = "spearman") {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method = method, use = "pairwise.complete.obs")
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  list(p = p.mat)
}

p_values_corr <- cor.mtest(corr_df, method = 'spearman')

# Save correlation plot
pdf(paste0(base_path, "SuppFig3_PET_Correlation_Matrix.pdf"), width = 10, height = 10)
corrplot(corr_matrix, 
         method = "color",
         type = "upper",
         addCoef.col = "black",
         number.cex = 0.8,
         p.mat = p_values_corr$p,
         sig.level = 0.05,
         insig = "blank",
         tl.cex = 1,
         tl.col = 'black',
         tl.srt = 45,
         col = colorRampPalette(c("#4575b4", "white", "#d73027"))(200),
         title = "Supplemental Figure 3. PET Metrics Correlation Matrix",
         mar = c(0, 0, 2, 0))
dev.off()

png(paste0(base_path, "SuppFig3_PET_Correlation_Matrix.png"), width = 10, height = 10, units = "in", res = 300)
corrplot(corr_matrix, 
         method = "color",
         type = "upper",
         addCoef.col = "black",
         number.cex = 0.8,
         p.mat = p_values_corr$p,
         sig.level = 0.05,
         insig = "blank",
         tl.cex = 1,
         tl.col = 'black',
         tl.srt = 45,
         col = colorRampPalette(c("#4575b4", "white", "#d73027"))(200),
         title = "Supplemental Figure 3. PET Metrics Correlation Matrix",
         mar = c(0, 0, 2, 0))
dev.off()

print("Supplemental Figure 3 saved!")

########################################################################
# SUPPLEMENTAL FIGURE 4: Arteriolohyalinosis vs PET Metrics
########################################################################

dat_clinical <- fread("C:/Users/netio/Downloads/UACR_Allparticipants_forGBM.csv")
dat_clinical <- dat_clinical[-str_which(dat_clinical$record_id, '-O')]

combined_df_clinical <- dat_clinical %>% 
  dplyr::select(record_id, avg_c_k2, avg_c_f, avg_c_k2_f, arteriolohyalinosis) %>%
  filter(!is.na(arteriolohyalinosis) & arteriolohyalinosis != '') %>%
  mutate(arteriolohyalinosis = factor(arteriolohyalinosis, 
                                      levels = c("no", "mild", "severe")))

PET_traits <- c('Cortical K2', 'Cortical F', 'Cortical K2/F')
PET_traits_small <- c('avg_c_k2', 'avg_c_f', 'avg_c_k2_f')
colors_arterio <- c("no" = "#00008B", "mild" = "#FFA500", "severe" = "#8B0000")

arterio_hyalin_plots <- list()
for(i in 1:length(PET_traits)){
  tmp_df <- combined_df_clinical
  iter <- which(names(tmp_df) == PET_traits_small[i])
  names(tmp_df)[iter] <- 'Variable'
  
  plot_data <- tmp_df %>% filter(!is.na(Variable))
  
  if(nrow(plot_data) > 0 && length(unique(plot_data$arteriolohyalinosis)) > 1){
    arterio_hyalin_plots[[i]] <- ggplot(plot_data, 
                                        aes(x = arteriolohyalinosis, y = Variable, 
                                            fill = arteriolohyalinosis)) +
      geom_boxplot() +
      scale_fill_manual(values = colors_arterio) +
      stat_compare_means(method = "kruskal.test", label = "p.format", size = 5) +
      theme_classic(base_size = 18) +
      theme(
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        plot.tag = element_text(size = 16, face = "bold")
      ) +
      labs(x = 'Arteriolohyalinosis Severity', y = PET_traits[i], 
           fill = 'Severity', tag = LETTERS[i])
  }
}

# Only combine plots that exist
valid_plots <- arterio_hyalin_plots[!sapply(arterio_hyalin_plots, is.null)]

if(length(valid_plots) == 3){
  p_supp4 <- (valid_plots[[1]] | valid_plots[[2]] | valid_plots[[3]]) +
    plot_annotation(
      title = "Supplemental Figure 4.",
      subtitle = 'Arteriolohyalinosis Severity and PET Imaging Biomarkers',
      theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
                    plot.subtitle = element_text(size = 14, hjust = 0))
    )
  
  ggsave(paste0(base_path, "SuppFig4_Arteriolohyalinosis_PET.pdf"), 
         plot = p_supp4, width = 18, height = 8)
  ggsave(paste0(base_path, "SuppFig4_Arteriolohyalinosis_PET.png"), 
         plot = p_supp4, width = 18, height = 8, dpi = 300)
  ggsave(paste0(base_path, "SuppFig4_Arteriolohyalinosis_PET.tiff"), 
         plot = p_supp4, width = 18, height = 8, dpi = 300, compression = "lzw")
  
  print("Supplemental Figure 4 saved!")
} else {
  print("Warning: Not enough data for Supplemental Figure 4")
}

########################################################################
# SUPPLEMENTAL FIGURE 5: Gene-level Statistics Tables for TCA/OxPhos
########################################################################

lc_files <- list.files('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', pattern='csv')
celltypes <- c('PT', 'PT-S1/S2', 'PT-S3', 'aPT')

# Collect all gene statistics
all_gene_stats <- data.frame()

for(i in 1:length(celltypes)){
  celltype <- celltypes[i]
  celltype2 <- str_replace_all(str_replace_all(celltype, "/", "_"), "-", "_")
  
  tmp_lc <- lc_files[str_which(lc_files, pattern = paste0('cycle_', celltype2, '_cells'))]
  
  # TCA genes
  tmp_lc_tca <- tmp_lc[str_which(tmp_lc, pattern = 'TCA')]
  if(length(tmp_lc_tca) > 0){
    tca_lc <- fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', tmp_lc_tca))
    tca_stats <- tca_lc %>% 
      dplyr::select(gene, 
                    logFC = logFC_groupType_2_Diabetes,
                    SE = any_of(c("se_groupType_2_Diabetes", "SE")),
                    pvalue = any_of(c("p_groupType_2_Diabetes", "pvalue"))) %>%
      mutate(pathway = "TCA Cycle",
             cell_type = celltype,
             significant = ifelse(pvalue < 0.05, "Yes", "No"))
    all_gene_stats <- bind_rows(all_gene_stats, tca_stats)
  }
  
  # OxPhos genes
  tmp_lc_oxphos <- tmp_lc[str_which(tmp_lc, pattern = 'PHOS_')]
  if(length(tmp_lc_oxphos) > 0){
    oxphos_lc <- fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', tmp_lc_oxphos))
    oxphos_stats <- oxphos_lc %>% 
      dplyr::select(gene, 
                    logFC = logFC_groupType_2_Diabetes,
                    SE = any_of(c("se_groupType_2_Diabetes", "SE")),
                    pvalue = any_of(c("p_groupType_2_Diabetes", "pvalue"))) %>%
      mutate(pathway = "OxPhos",
             cell_type = celltype,
             significant = ifelse(pvalue < 0.05, "Yes", "No"))
    all_gene_stats <- bind_rows(all_gene_stats, oxphos_stats)
  }
}

# Save comprehensive gene statistics table
write.csv(all_gene_stats, 
          paste0(base_path, 'SuppTable_Gene_Statistics_All.csv'), 
          row.names = FALSE)

# Create summary by pathway and cell type
gene_summary <- all_gene_stats %>%
  group_by(pathway, cell_type) %>%
  summarise(
    n_genes = n(),
    n_significant = sum(significant == "Yes"),
    pct_significant = round(100 * n_significant / n_genes, 1),
    n_upregulated = sum(logFC > 0 & significant == "Yes"),
    n_downregulated = sum(logFC < 0 & significant == "Yes"),
    mean_logFC = round(mean(logFC, na.rm = TRUE), 4),
    .groups = 'drop'
  )

write.csv(gene_summary, 
          paste0(base_path, 'SuppTable_Gene_Summary.csv'), 
          row.names = FALSE)

print("Supplemental Tables saved!")

cat("\n=== Gene Expression Summary ===\n")
print(gene_summary)

########################################################################
# SUPPLEMENTAL FIGURE 6: Demographics Table for All Cohorts
########################################################################

# Demographics for PET analysis (all 3 groups)
combined_df_demo <- dat_all %>%
  mutate(
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    hba1c = as.numeric(hba1c),
    acr_u = as.numeric(acr_u),
    egfr = as.numeric(egfr),
    sex = as.factor(sex),
    race_ethnicity = as.factor(race_ethnicity),
    study = as.factor(study),
    group = factor(group, levels = c("Lean Control", "Obese Control", "Type 2 Diabetes"))
  )

desc_table_3groups <- combined_df_demo %>%
  dplyr::select(age, sex, race_ethnicity, bmi, hba1c, study, group, acr_u, egfr) %>%
  tbl_summary(
    by = group,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      hba1c ~ "continuous",
      acr_u ~ "continuous",
      egfr ~ "continuous",
      sex ~ "categorical",
      race_ethnicity ~ "categorical",
      study ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      acr_u ~ 1,
      egfr ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      sex ~ "Sex", 
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      study ~ "Study",
      acr_u ~ "UACR, mg/g",
      egfr ~ "eGFR, ml/min/1.73m²"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "kruskal.test",
    all_categorical() ~ "fisher.test"
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Group**")

desc_table_3groups %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave(paste0(base_path, "SuppTable_Demographics_3Groups.png"), 
         vwidth = 1400, vheight = 900)

print("Supplemental demographics table saved!")

########################################################################
# COMBINE ALL SUPPLEMENTAL PDFs
########################################################################

library(qpdf)

supp_pdf_files <- c(
  paste0(base_path, "SuppFig1_PET_3Groups.pdf"),
  paste0(base_path, "SuppFig2_Medullary_PET.pdf"),
  paste0(base_path, "SuppFig3_PET_Correlation_Matrix.pdf"),
  paste0(base_path, "SuppFig4_Arteriolohyalinosis_PET.pdf")
)

existing_supp_files <- supp_pdf_files[file.exists(supp_pdf_files)]

if(length(existing_supp_files) > 0) {
  pdf_combine(input = existing_supp_files, 
              output = paste0(base_path, "Combined_Supplemental_Figures.pdf"))
  print("Combined Supplemental PDF saved!")
} else {
  print("No supplemental PDF files found to combine")
}

print("All supplemental figures complete!")

