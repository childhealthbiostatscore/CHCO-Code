########################################################################
# ROCKIES Publication Figures - Updated Organization
# Figure 2: PET Comparisons (T2D vs HC)
# Figure 3: UACR Correlation
# Figure 4: GBM PET Comparisons
# Figure 5: Arteriosclerosis PET Comparisons
# Figure 6: scRNAseq TCA Cycle
# Figure 7: scRNAseq OxPhos
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
# DATA PREPARATION
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

# Remove existing PET average columns if they exist to avoid duplicates
dat <- dat %>% dplyr::select(-any_of(c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 'avg_c_k2_f', 'avg_m_k2_f')))

tmp_results <- PET_avg(dat)
dat_results <- dat %>% bind_cols(tmp_results)
dat_results <- dat_results %>% filter(!is.na(avg_c_k2))
dat_results <- dat_results %>% filter(group %in% c('Lean Control', 'Type 2 Diabetes'))

# Get SGLT2i information
RH <- fread('C:/Users/netio/Documents/UofW/Rockies/RENALHEIR-SGLT2.csv')
names(RH) <- c('Subject', 'rep_instr', 'rep_inst', 'SGLT2')
RH2 <- fread('C:/Users/netio/Documents/UofW/Rockies/RenalHEIRitage-SGLT2Use.csv')
names(RH2) <- c('Subject', 'event', 'rep_instr', 'rep_inst', 'mrn', 'SGLT2', 'SGLT2_ever')
RH2 <- RH2 %>% filter(!is.na(mrn))
improve <- fread('C:/Users/netio/Downloads/IMPROVET2D-SGLT2i_DATA_LABELS_2025-08-25_0938.csv')
names(improve)[5] <- 'SGLT2'
names(improve)[1] <- 'record_id'
improve <- improve %>% filter(!is.na(SGLT2)) %>% filter(SGLT2 != '')

dat2 <- dat_results
dat2$group2 <- NA
need_med_info <- dat2 %>% filter(is.na(group2))
improve_small <- improve %>% filter(record_id %in% need_med_info$record_id)
RH_small <- RH %>% filter(Subject %in% need_med_info$record_id)
RH2_small <- RH2 %>% filter(mrn %in% need_med_info$mrn)

for(i in c(1:nrow(RH_small))){
  if(nrow(RH_small) == 0) next
  if(RH_small$SGLT2[i] == 'No'){
    dat2$group2[which(dat2$record_id == RH_small$Subject[i])] <- 'T2D-No SGLTi2'
    dat2$epic_sglti2_1[which(dat2$record_id == RH_small$Subject[i])] <- 'No'
  }else if(RH_small$SGLT2[i] == 'Yes'){
    dat2$group2[which(dat2$record_id == RH_small$Subject[i])] <- 'T2D-SGLTi2'
    dat2$epic_sglti2_1[which(dat2$record_id == RH_small$Subject[i])] <- 'Yes'
  }
}

for(i in c(1:nrow(RH2_small))){
  if(nrow(RH2_small) == 0) next
  if(RH2_small$SGLT2[i] == 'No'){
    dat2$group2[which(dat2$mrn == RH2_small$mrn[i])] <- 'T2D-No SGLTi2'
    dat2$epic_sglti2_1[which(dat2$mrn == RH2_small$mrn[i])] <- 'No'
  }else if(RH2_small$SGLT2[i] == 'Yes'){
    dat2$group2[which(dat2$mrn == RH2_small$mrn[i])] <- 'T2D-SGLTi2'
    dat2$epic_sglti2_1[which(dat2$mrn == RH2_small$mrn[i])] <- 'Yes'
  }
}

for(i in c(1:nrow(improve_small))){
  if(nrow(improve_small) == 0) next
  if(improve_small$SGLT2[i] == 'No'){
    dat2$group2[which(dat2$record_id == improve_small$record_id[i])] <- 'T2D-No SGLTi2'
    dat2$epic_sglti2_1[which(dat2$record_id == improve_small$record_id[i])] <- 'No'
  }else if(improve_small$SGLT2[i] == 'Yes'){
    dat2$group2[which(dat2$record_id == improve_small$record_id[i])] <- 'T2D-SGLTi2'
    dat2$epic_sglti2_1[which(dat2$record_id == improve_small$record_id[i])] <- 'Yes'
  }
}

dat2$epic_sglti2_1[which(dat2$group == 'Lean Control')] <- 'No'
dat2 <- dat2 %>% filter(epic_sglti2_1 != 'Yes')

########################################################################
# FIGURE 2: PET Comparisons (T2D vs HC)
########################################################################

aim1_df <- dat2 %>% dplyr::select(record_id, group, avg_c_f, avg_c_k2, avg_c_k2_f)

aim1_long <- aim1_df %>%
  pivot_longer(cols = c(avg_c_f, avg_c_k2, avg_c_k2_f), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("avg_c_f", "avg_c_k2", "avg_c_k2_f"),
                         labels = c("Cortical F", "Cortical K2", "Cortical K2/F")))

stat_test <- aim1_long %>% group_by(metric) %>%
  pairwise_wilcox_test(value ~ group, p.adjust.method = "none") %>%
  add_significance() %>% add_xy_position(x = "metric", dodge = 0.8)

stat_test_adjusted <- stat_test %>% mutate(y.position = c(2.9, 0.25, 0.28))

p_fig2 <- ggplot(aim1_long, aes(x = metric, y = value)) +
  geom_boxplot(aes(fill = group), color = "black", alpha = 0.8, outlier.shape = NA,
               linewidth = 0.5, position = position_dodge(width = 0.8)) +
  geom_point(aes(fill = group), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, size = 3, shape = 21, color = "black", stroke = 0.3) +
  scale_fill_manual(values = c("#c2dfe3", "#fcb1a6"), labels = c("Lean Control", "Type 2 Diabetes")) +
  labs(x = "PET Metrics", y = "Value", fill = "Group",
       title = "Figure 2.", subtitle = "Global PET Imaging Metrics by Group") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title = element_text(size = 22, face = "bold", color = "black"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0, color = "black"),
        plot.subtitle = element_text(size = 16, hjust = 0, color = "black", margin = margin(b = 10)),
        legend.position = "bottom", legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 18, face = "bold", color = "black"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

p_fig2_stats <- p_fig2 + stat_pvalue_manual(stat_test_adjusted, label = "p.adj.signif", 
                                            tip.length = 0.01, hide.ns = TRUE, size = 6)
p_fig2_final <- p_fig2_stats + scale_y_break(c(0.3, 0.70), scales = 2)

ggsave(paste0(base_path, 'Figure2_PET_Comparisons.pdf'), plot = p_fig2_final, width = 12, height = 14, device = cairo_pdf, dpi = 300)
ggsave(paste0(base_path, 'Figure2_PET_Comparisons.png'), plot = p_fig2_final, width = 12, height = 14, dpi = 300)
ggsave(paste0(base_path, 'Figure2_PET_Comparisons.tiff'), plot = p_fig2_final, width = 12, height = 14, dpi = 300, compression = "lzw")
print("Figure 2 saved!")

########################################################################
# FIGURE 3: UACR Correlation
########################################################################

# Need to recalculate PET averages for the full dataset (not just Lean/T2D)
dat_for_uacr <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))

dat_for_uacr <- dat_for_uacr %>% dplyr::select(-any_of(c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 'avg_c_k2_f', 'avg_m_k2_f')))
tmp_results_uacr <- PET_avg(dat_for_uacr)
dat_for_uacr <- dat_for_uacr %>% bind_cols(tmp_results_uacr)

dat_uacr <- dat_for_uacr %>% filter(!is.na(avg_c_k2)) %>% 
  filter(group %in% c('Lean Control', 'Obese Control', 'Type 2 Diabetes')) %>% 
  mutate(avg_c_k2_f = avg_c_k2 / avg_c_f)
dat_uacr <- dat_uacr %>% filter(record_id != 'CRC-55')
dat_uacr$group[which(dat_uacr$record_id == 'RH2-39-O')] <- 'Obese Control'

combined_df_uacr <- dat_uacr %>% dplyr::select(avg_c_k2, avg_c_f, avg_c_k2_f, acr_u)
combined_df_corr <- cor(combined_df_uacr, use = 'pairwise.complete.obs', method = 'spearman')

cor.mtest <- function(mat, method = "spearman") {
  mat <- as.matrix(mat); n <- ncol(mat); p.mat <- matrix(NA, n, n); diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method = method, use = "pairwise.complete.obs")
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  list(p = p.mat)
}

p_values <- cor.mtest(combined_df_uacr, method = 'spearman')
corr_subset <- as.matrix(combined_df_corr[4, 1:3, drop = FALSE])
p_subset <- as.matrix(p_values$p[4, 1:3, drop = FALSE])
rownames(corr_subset) <- c('UACR'); colnames(corr_subset) <- c('Cortical K2', 'Cortical F', 'Cortical K2/F')
rownames(p_subset) <- rownames(corr_subset); colnames(p_subset) <- colnames(corr_subset)

pdf(paste0(base_path, "Figure3_UACR_Correlation.pdf"), width = 12, height = 6)
par(mar = c(2, 2, 4, 2))
corrplot(corr_subset, method = "color", p.mat = p_subset, sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 2, insig = "label_sig", number.cex = 1.5, tl.cex = 1.5, tl.col = 'black',
         cl.cex = 1.2, mar = c(0, 0, 2, 0), title = "Figure 3. UACR Correlation with PET Metrics")
dev.off()

png(paste0(base_path, "Figure3_UACR_Correlation.png"), width = 12, height = 6, units = "in", res = 300)
par(mar = c(2, 2, 4, 2))
corrplot(corr_subset, method = "color", p.mat = p_subset, sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 2, insig = "label_sig", number.cex = 1.5, tl.cex = 1.5, tl.col = 'black',
         cl.cex = 1.2, mar = c(0, 0, 2, 0), title = "Figure 3. UACR Correlation with PET Metrics")
dev.off()

tiff(paste0(base_path, "Figure3_UACR_Correlation.tiff"), width = 12, height = 6, units = "in", res = 300, compression = "lzw")
par(mar = c(2, 2, 4, 2))
corrplot(corr_subset, method = "color", p.mat = p_subset, sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 2, insig = "label_sig", number.cex = 1.5, tl.cex = 1.5, tl.col = 'black',
         cl.cex = 1.2, mar = c(0, 0, 2, 0), title = "Figure 3. UACR Correlation with PET Metrics")
dev.off()
print("Figure 3 saved!")

########################################################################
# FIGURE 4: GBM PET Comparisons
########################################################################

dat_clinical <- fread("C:/Users/netio/Downloads/UACR_Allparticipants_forGBM.csv")
dat_clinical <- dat_clinical[-str_which(dat_clinical$record_id, '-O')]
combined_df_clinical <- dat_clinical %>% 
  dplyr::select(record_id, avg_c_k2, avg_c_f, avg_c_k2_f, arteriosclerosis, `GBM thickness`) %>% 
  mutate(GBM_thickness = `GBM thickness`)

PET_traits <- c('Cortical K2', 'Cortical F', 'Cortical K2/F')
PET_traits_small <- c('avg_c_k2', 'avg_c_f', 'avg_c_k2_f')
colors <- c("no" = "#00008B", "yes" = "#8B0000")

gbm_plots <- list()
for(i in 1:length(PET_traits)){
  tmp_df <- combined_df_clinical; iter <- which(names(tmp_df) == PET_traits_small[i])
  names(tmp_df)[iter] <- 'Variable'
  plot_data <- tmp_df %>% filter(GBM_thickness != '' & !is.na(GBM_thickness))
  gbm_plots[[i]] <- ggplot(plot_data, aes(x = as.character(GBM_thickness), y = Variable, fill = as.character(GBM_thickness))) +
    geom_boxplot() + scale_fill_manual(values = colors) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 6) +
    theme_classic(base_size = 20) +
    theme(axis.title = element_text(size = 22, face = "bold"), axis.text = element_text(size = 20),
          legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 18),
          plot.tag = element_text(size = 16, face = "bold")) +
    labs(x = 'GBM Thickening', y = PET_traits[i], fill = 'GBM Thickening', tag = LETTERS[i])
}

p_fig4 <- (gbm_plots[[1]] | gbm_plots[[2]] | gbm_plots[[3]]) +
  plot_annotation(title = "Figure 4.", subtitle = 'GBM Thickening and PET Imaging Biomarkers',
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
                                plot.subtitle = element_text(size = 16, hjust = 0)))

ggsave(paste0(base_path, "Figure4_GBM_PET.pdf"), plot = p_fig4, width = 18, height = 8)
ggsave(paste0(base_path, "Figure4_GBM_PET.png"), plot = p_fig4, width = 18, height = 8, dpi = 300)
ggsave(paste0(base_path, "Figure4_GBM_PET.tiff"), plot = p_fig4, width = 18, height = 8, dpi = 300, compression = "lzw")
print("Figure 4 saved!")

########################################################################
# FIGURE 5: Arteriosclerosis PET Comparisons
########################################################################

arterio_plots <- list()
for(i in 1:length(PET_traits)){
  tmp_df <- combined_df_clinical; iter <- which(names(tmp_df) == PET_traits_small[i])
  names(tmp_df)[iter] <- 'Variable'
  plot_data <- tmp_df %>% filter(arteriosclerosis != '' & !is.na(arteriosclerosis))
  arterio_plots[[i]] <- ggplot(plot_data, aes(x = as.character(arteriosclerosis), y = Variable, fill = as.character(arteriosclerosis))) +
    geom_boxplot() + scale_fill_manual(values = colors) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 6) +
    theme_classic(base_size = 20) +
    theme(axis.title = element_text(size = 22, face = "bold"), axis.text = element_text(size = 20),
          legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 18),
          plot.tag = element_text(size = 16, face = "bold")) +
    labs(x = 'Arteriosclerosis', y = PET_traits[i], fill = 'Arteriosclerosis', tag = LETTERS[i])
}

p_fig5 <- (arterio_plots[[1]] | arterio_plots[[2]] | arterio_plots[[3]]) +
  plot_annotation(title = "Figure 5.", subtitle = 'Arteriosclerosis and PET Imaging Biomarkers',
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
                                plot.subtitle = element_text(size = 16, hjust = 0)))

ggsave(paste0(base_path, "Figure5_Arteriosclerosis_PET.pdf"), plot = p_fig5, width = 18, height = 8)
ggsave(paste0(base_path, "Figure5_Arteriosclerosis_PET.png"), plot = p_fig5, width = 18, height = 8, dpi = 300)
ggsave(paste0(base_path, "Figure5_Arteriosclerosis_PET.tiff"), plot = p_fig5, width = 18, height = 8, dpi = 300, compression = "lzw")
print("Figure 5 saved!")

########################################################################
# FIGURE 6: scRNAseq TCA Cycle
########################################################################

lc_files <- list.files('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', pattern='csv')
celltypes <- c('PT', 'PT-S1/S2', 'PT-S3', 'aPT')
tca_plots <- list()

for(i in 1:length(celltypes)){
  celltype <- celltypes[i]
  if(celltype == 'PT'){ legend_position <- 'bottom'; axis_angle <- 0; hjust_val <- 0.5; asterisk_size <- 8
  } else { legend_position <- 'none'; axis_angle <- 45; hjust_val <- 1; asterisk_size <- 6 }
  
  celltype2 <- str_replace_all(str_replace_all(celltype, "/", "_"), "-", "_")
  tmp_lc <- lc_files[str_which(lc_files, pattern = paste0('cycle_', celltype2, '_cells'))]
  tmp_lc_tca <- tmp_lc[str_which(tmp_lc, pattern = 'TCA')]
  tca_lc <- fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', tmp_lc_tca))
  tca_lc <- tca_lc %>% dplyr::select(gene, logFC_lc = logFC_groupType_2_Diabetes, 
                                     pvalue_lc = any_of(c("p_groupType_2_Diabetes", "pvalue")))
  tca_full_sig <- tca_lc %>% mutate(direction = ifelse(logFC_lc < 0, 'Down', ifelse(logFC_lc > 0, 'Up', NA)))
  
  if(nrow(tca_full_sig) > 0){
    plot_df <- tca_full_sig %>% gather(key = 'metric_condition', value = 'value', -gene, -direction) %>% 
      separate(metric_condition, into = c('metric', 'condition'), sep = '_') %>% 
      spread(key = metric, value = value) %>% mutate(Significance = ifelse(pvalue < 0.05, '*', '')) %>% 
      filter(condition == 'lc')
    
    tca_plots[[celltype]] <- ggplot(plot_df, aes(x=gene, y=logFC, fill=direction)) +
      geom_bar(stat='identity', position='dodge', color='black', linewidth=0.3) +
      geom_text(aes(label = Significance), position = position_dodge(width = 0.9), vjust = -0.5, size = asterisk_size) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
      labs(x = NULL, y = 'Log2 Fold Change', subtitle = celltype) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      scale_fill_manual(name = 'Regulation', labels = c('Up' = 'Up-regulated', 'Down' = 'Down-regulated'),
                        values = c('Down' = '#d73027', 'Up' = '#4575b4')) +
      theme_classic(base_size = 11) +
      theme(axis.text.x = element_text(angle = axis_angle, hjust = hjust_val, size = 10, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title = element_text(size = 11, face = "bold", color = "black"),
            legend.position = legend_position,
            plot.subtitle = element_text(size = 11, hjust = 0.5, face = "bold"),
            panel.background = element_rect(fill = "white", color = NA))
  }
}

legend_tca <- get_legend(tca_plots[['PT']])
tca_plots[['PT']] <- tca_plots[['PT']] + theme(legend.position = "none") + labs(tag = "A") + theme(plot.tag = element_text(size = 14, face = "bold"))
tca_plots[['PT-S1/S2']] <- tca_plots[['PT-S1/S2']] + labs(tag = "B") + theme(plot.tag = element_text(size = 14, face = "bold"))
tca_plots[['PT-S3']] <- tca_plots[['PT-S3']] + labs(tag = "C") + theme(plot.tag = element_text(size = 14, face = "bold"))
tca_plots[['aPT']] <- tca_plots[['aPT']] + labs(tag = "D") + theme(plot.tag = element_text(size = 14, face = "bold"))

tca_combined <- (tca_plots[['PT']]) / (tca_plots[['PT-S1/S2']] | tca_plots[['PT-S3']] | tca_plots[['aPT']]) +
  plot_layout(heights = c(1.3, 1.3)) +
  plot_annotation(title = 'Figure 6.', subtitle = 'TCA cycle gene expression in proximal tubule cells',
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
                                plot.subtitle = element_text(size = 16, hjust = 0, margin = margin(b = 10))))

tca_final <- tca_combined + inset_element(legend_tca, left = 0.43, bottom = -0.02, right = 0.57, top = 0.02, align_to = 'full')

ggsave(paste0(base_path, 'Figure6_TCA_Cycle.pdf'), plot = tca_final, width = 16, height = 12, device = cairo_pdf, dpi = 300)
ggsave(paste0(base_path, 'Figure6_TCA_Cycle.png'), plot = tca_final, width = 16, height = 12, dpi = 300)
ggsave(paste0(base_path, 'Figure6_TCA_Cycle.tiff'), plot = tca_final, width = 16, height = 12, dpi = 300, compression = "lzw")
print("Figure 6 saved!")

########################################################################
# FIGURE 7: scRNAseq OxPhos
########################################################################

oxphos_plots <- list()

for(i in 1:length(celltypes)){
  celltype <- celltypes[i]
  if(celltype == 'PT'){ legend_position <- 'bottom'; axis_angle <- 0; hjust_val <- 0.5; asterisk_size <- 8
  } else { legend_position <- 'none'; axis_angle <- 45; hjust_val <- 1; asterisk_size <- 6 }
  
  celltype2 <- str_replace_all(str_replace_all(celltype, "/", "_"), "-", "_")
  tmp_lc <- lc_files[str_which(lc_files, pattern = paste0('cycle_', celltype2, '_cells'))]
  tmp_lc_oxphos <- tmp_lc[str_which(tmp_lc, pattern = 'PHOS_')]
  oxphos_lc <- fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', tmp_lc_oxphos))
  oxphos_lc <- oxphos_lc %>% dplyr::select(gene, logFC_lc = logFC_groupType_2_Diabetes, 
                                           pvalue_lc = any_of(c("p_groupType_2_Diabetes", "pvalue")))
  oxphos_full_sig <- oxphos_lc %>% mutate(direction = ifelse(logFC_lc < 0, 'Down', ifelse(logFC_lc > 0, 'Up', NA)))
  
  if(nrow(oxphos_full_sig) > 0){
    plot_df <- oxphos_full_sig %>% gather(key = 'metric_condition', value = 'value', -gene, -direction) %>% 
      separate(metric_condition, into = c('metric', 'condition'), sep = '_') %>% 
      spread(key = metric, value = value) %>% mutate(Significance = ifelse(pvalue < 0.05, '*', '')) %>% 
      filter(condition == 'lc')
    
    oxphos_plots[[celltype]] <- ggplot(plot_df, aes(x=gene, y=logFC, fill=direction)) +
      geom_bar(stat='identity', position='dodge', color='black', linewidth=0.3) +
      geom_text(aes(label = Significance), position = position_dodge(width = 0.9), vjust = -0.5, size = asterisk_size) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
      labs(x = NULL, y = 'Log2 Fold Change', subtitle = celltype) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      scale_fill_manual(name = 'Regulation', labels = c('Up' = 'Up-regulated', 'Down' = 'Down-regulated'),
                        values = c('Down' = '#d73027', 'Up' = '#4575b4')) +
      theme_classic(base_size = 11) +
      theme(axis.text.x = element_text(angle = axis_angle, hjust = hjust_val, size = 10, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title = element_text(size = 11, face = "bold", color = "black"),
            legend.position = legend_position,
            plot.subtitle = element_text(size = 11, hjust = 0.5, face = "bold"),
            panel.background = element_rect(fill = "white", color = NA))
  }
}

legend_oxphos <- get_legend(oxphos_plots[['PT']])
oxphos_plots[['PT']] <- oxphos_plots[['PT']] + theme(legend.position = "none") + labs(tag = "A") + theme(plot.tag = element_text(size = 14, face = "bold"))
oxphos_plots[['PT-S1/S2']] <- oxphos_plots[['PT-S1/S2']] + labs(tag = "B") + theme(plot.tag = element_text(size = 14, face = "bold"))
oxphos_plots[['PT-S3']] <- oxphos_plots[['PT-S3']] + labs(tag = "C") + theme(plot.tag = element_text(size = 14, face = "bold"))
oxphos_plots[['aPT']] <- oxphos_plots[['aPT']] + labs(tag = "D") + theme(plot.tag = element_text(size = 14, face = "bold"))

oxphos_combined <- (oxphos_plots[['PT']]) / (oxphos_plots[['PT-S1/S2']] | oxphos_plots[['PT-S3']] | oxphos_plots[['aPT']]) +
  plot_layout(heights = c(1.3, 1.3)) +
  plot_annotation(title = 'Figure 7.', subtitle = 'Oxidative phosphorylation gene expression in proximal tubule cells',
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
                                plot.subtitle = element_text(size = 16, hjust = 0, margin = margin(b = 10))))

oxphos_final <- oxphos_combined + inset_element(legend_oxphos, left = 0.43, bottom = -0.02, right = 0.57, top = 0.02, align_to = 'full')

ggsave(paste0(base_path, 'Figure7_OxPhos.pdf'), plot = oxphos_final, width = 16, height = 12, device = cairo_pdf, dpi = 300)
ggsave(paste0(base_path, 'Figure7_OxPhos.png'), plot = oxphos_final, width = 16, height = 12, dpi = 300)
ggsave(paste0(base_path, 'Figure7_OxPhos.tiff'), plot = oxphos_final, width = 16, height = 12, dpi = 300, compression = "lzw")
print("Figure 7 saved!")

########################################################################
# COMBINE ALL PDFs
########################################################################

library(qpdf)

pdf_files <- c(
  paste0(base_path, "Figure2_PET_Comparisons.pdf"),
  paste0(base_path, "Figure3_UACR_Correlation.pdf"),
  paste0(base_path, "Figure4_GBM_PET.pdf"),
  paste0(base_path, "Figure5_Arteriosclerosis_PET.pdf"),
  paste0(base_path, "Figure6_TCA_Cycle.pdf"),
  paste0(base_path, "Figure7_OxPhos.pdf")
)

existing_files <- pdf_files[file.exists(pdf_files)]

if(length(existing_files) > 0) {
  pdf_combine(input = existing_files, output = paste0(base_path, "Combined_All_Figures.pdf"))
  print("Combined PDF saved!")
} else {
  print("No PDF files found to combine")
}

print("All figures complete!")





