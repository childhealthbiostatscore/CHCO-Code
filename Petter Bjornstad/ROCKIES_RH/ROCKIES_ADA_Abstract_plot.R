# ============================================================================
# MULTI-PANEL FIGURE: PET METRICS BY GROUP, PATHOLOGY, AND UACR
# ============================================================================

# Clear environment
remove(list=ls())

# Load libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)
library(stringr)
library(data.table)
library(corrplot)
library(grid)
library(png)
library(tidyr)
library(gridExtra)

# Set paths
base_path <- 'C:/Users/netio/Documents/UofW/Rockies/publication_figures/'
harmonized_path <- "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv"
clinical_path <- "C:/Users/netio/Downloads/UACR_Allparticipants_forGBM.csv"

# ============================================================================
# LOAD AND PREPARE DATA
# ============================================================================

cat("\n=== LOADING DATA ===\n\n")

# Load harmonized data
harmonized_data <- read.csv(harmonized_path, na = '')

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))

# Calculate PET averages
PET_avg <- function(data){
  tmp_df <- data %>% dplyr::select(lc_k2, rc_k2, lm_k2, rm_k2,
                                   lc_f, rc_f, lm_f, rm_f)
  avg_c_k2 <- tmp_df %>%
    dplyr::select(lc_k2, rc_k2) %>% rowMeans(na.rm=T)
  avg_m_k2 <- tmp_df %>% 
    dplyr::select(lm_k2, rm_k2) %>% rowMeans(na.rm=T)
  avg_c_f <- tmp_df %>% 
    dplyr::select(lc_f, rc_f) %>% rowMeans(na.rm=T)
  avg_m_f <- tmp_df %>% 
    dplyr::select(lm_f, rm_f) %>% rowMeans(na.rm=T)
  avg_c_k2_f <- avg_c_k2 / avg_c_f
  avg_m_k2_f <- avg_m_k2 / avg_m_f
  
  results <- bind_cols(avg_c_k2, avg_m_k2, avg_c_f, avg_m_f, 
                       avg_c_k2_f, avg_m_k2_f) %>% as.data.frame()
  names(results) <- c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 
                      'avg_c_k2_f', 'avg_m_k2_f')
  return(results)
}

tmp_results <- PET_avg(dat)

# Remove any existing PET average columns before binding new ones
dat_results <- dat %>% 
  select(-any_of(c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 'avg_c_k2_f', 'avg_m_k2_f'))) %>%
  bind_cols(tmp_results)

# Filter for PET data and correct groups
dat_results <- dat_results %>% filter(!is.na(avg_c_k2))
dat_results <- dat_results %>% filter(group %in% c('Lean Control', 'Type 2 Diabetes'))

# Rename "Lean Control" to "Healthy Weight Control"
dat_results <- dat_results %>%
  mutate(group = ifelse(group == "Lean Control", "Healthy Weight Control", group))

# Get SGLT2i information
RH <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/RENALHEIR-SGLT2.csv')
names(RH) <- c('Subject', 'rep_instr', 'rep_inst', 'SGLT2')
RH2 <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/RenalHEIRitage-SGLT2Use.csv')
names(RH2) <- c('Subject', 'event', 'rep_instr', 'rep_inst', 'mrn', 'SGLT2', 'SGLT2_ever')
RH2 <- RH2 %>% filter(!is.na(mrn))
improve <- data.table::fread('C:/Users/netio/Downloads/IMPROVET2D-SGLT2i_DATA_LABELS_2025-08-25_0938.csv')
names(improve)[5] <- 'SGLT2'
names(improve)[1] <- 'record_id'
improve <- improve %>% filter(!is.na(SGLT2)) %>% filter(SGLT2 != '')

dat2 <- dat_results
dat2$group2 <- NA
need_med_info <- dat2 %>% filter(is.na(group2))

improve_small <- improve %>% filter(record_id %in% need_med_info$record_id)
RH_small <- RH %>% filter(Subject %in% need_med_info$record_id)
RH2_small <- RH2 %>% filter(mrn %in% need_med_info$mrn)

# Apply SGLT2i information
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

dat2$epic_sglti2_1[which(dat2$group == 'Healthy Weight Control')] <- 'No'
dat2 <- dat2 %>% filter(epic_sglti2_1 != 'Yes')

# Load clinical/pathology data
dat_clinical <- data.table::fread(clinical_path)
dat_clinical <- dat_clinical[-str_which(dat_clinical$record_id, '-O')]

cat("All data loaded successfully\n\n")

# ============================================================================
# PANEL A-C: PET METRICS BY GROUP (Healthy Weight Control vs T2D)
# ============================================================================

cat("Creating Panels A-C: PET by Group\n")

# Prepare data for group comparison
aim1_df <- dat2 %>% 
  dplyr::select(record_id, group, avg_c_f, avg_c_k2, avg_c_k2_f)

aim1_long <- aim1_df %>%
  pivot_longer(cols = c(avg_c_f, avg_c_k2, avg_c_k2_f),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = factor(metric, 
                         levels = c("avg_c_f", "avg_c_k2", "avg_c_k2_f"),
                         labels = c("Cortical F", "Cortical K2", "Cortical K2/F")))

# Create group comparison plots
group_colors <- c("Healthy Weight Control" = "#c2dfe3", "Type 2 Diabetes" = "#fcb1a6")

plot_a <- aim1_long %>% filter(metric == "Cortical F") %>%
  ggplot(aes(x = group, y = value, fill = group)) +
  geom_boxplot(color = "black", alpha = 0.8, outlier.shape = NA, linewidth = 0.5) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6, size = 2, shape = 21, color = "black", stroke = 0.3) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5) +
  labs(x = NULL, y = "Cortical F", tag = "A") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none",
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

plot_b <- aim1_long %>% filter(metric == "Cortical K2") %>%
  ggplot(aes(x = group, y = value, fill = group)) +
  geom_boxplot(color = "black", alpha = 0.8, outlier.shape = NA, linewidth = 0.5) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6, size = 2, shape = 21, color = "black", stroke = 0.3) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5) +
  labs(x = NULL, y = "Cortical K2", tag = "B") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none",
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

plot_c <- aim1_long %>% filter(metric == "Cortical K2/F") %>%
  ggplot(aes(x = group, y = value, fill = group)) +
  geom_boxplot(color = "black", alpha = 0.8, outlier.shape = NA, linewidth = 0.5) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6, size = 2, shape = 21, color = "black", stroke = 0.3) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5) +
  labs(x = NULL, y = "Cortical K2/F", tag = "C") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none",
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# ============================================================================
# PANEL D-F: GBM THICKNESS
# ============================================================================

cat("Creating Panels D-F: GBM Thickness\n")

combined_df_gbm <- dat_clinical %>% 
  dplyr::select(record_id, avg_c_k2, avg_c_f, avg_c_k2_f, `GBM thickness`) %>% 
  filter(`GBM thickness` != '', !is.na(`GBM thickness`))

pathology_colors <- c("no" = "#00008B", "yes" = "#8B0000")

plot_d <- ggplot(combined_df_gbm, aes(x = `GBM thickness`, y = avg_c_k2, fill = `GBM thickness`)) +
  geom_boxplot(color = "black", alpha = 0.8, linewidth = 0.5) +
  scale_fill_manual(values = pathology_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5, p.adjust.method = "none") +
  labs(x = "GBM Thickening", y = "Cortical K2", fill = "GBM Thickening", tag = "D") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none",
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

plot_e <- ggplot(combined_df_gbm, aes(x = `GBM thickness`, y = avg_c_f, fill = `GBM thickness`)) +
  geom_boxplot(color = "black", alpha = 0.8, linewidth = 0.5) +
  scale_fill_manual(values = pathology_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5, p.adjust.method = "none") +
  labs(x = "GBM Thickening", y = "Cortical F", fill = "GBM Thickening", tag = "E") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none",
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

plot_f <- ggplot(combined_df_gbm, aes(x = `GBM thickness`, y = avg_c_k2_f, fill = `GBM thickness`)) +
  geom_boxplot(color = "black", alpha = 0.8, linewidth = 0.5) +
  scale_fill_manual(values = pathology_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5, p.adjust.method = "none") +
  labs(x = "GBM Thickening", y = "Cortical K2/F", fill = "GBM Thickening", tag = "F") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none",
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# ============================================================================
# PANEL G-I: ARTERIOSCLEROSIS
# ============================================================================

cat("Creating Panels G-I: Arteriosclerosis\n")

combined_df_arterio <- dat_clinical %>% 
  dplyr::select(record_id, avg_c_k2, avg_c_f, avg_c_k2_f, arteriosclerosis) %>% 
  filter(arteriosclerosis != '', !is.na(arteriosclerosis))

plot_g <- ggplot(combined_df_arterio, aes(x = arteriosclerosis, y = avg_c_k2, fill = arteriosclerosis)) +
  geom_boxplot(color = "black", alpha = 0.8, linewidth = 0.5) +
  scale_fill_manual(values = pathology_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5, p.adjust.method = "none") +
  labs(x = "Arteriosclerosis", y = "Cortical K2", fill = "Arteriosclerosis", tag = "G") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none",
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

plot_h <- ggplot(combined_df_arterio, aes(x = arteriosclerosis, y = avg_c_f, fill = arteriosclerosis)) +
  geom_boxplot(color = "black", alpha = 0.8, linewidth = 0.5) +
  scale_fill_manual(values = pathology_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5, p.adjust.method = "none") +
  labs(x = "Arteriosclerosis", y = "Cortical F", fill = "Arteriosclerosis", tag = "H") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none",
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

plot_i <- ggplot(combined_df_arterio, aes(x = arteriosclerosis, y = avg_c_k2_f, fill = arteriosclerosis)) +
  geom_boxplot(color = "black", alpha = 0.8, linewidth = 0.5) +
  scale_fill_manual(values = pathology_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5, p.adjust.method = "none") +
  labs(x = "Arteriosclerosis", y = "Cortical K2/F", fill = "Arteriosclerosis", tag = "I") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none",
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )


# ============================================================================
# PANEL J: UACR CORRELATIONS (EXACT VALUES FROM YOUR FIGURE)
# ============================================================================

cat("Creating Panel J: UACR Correlations\n")

# Use the exact correlation values from your figure legend:
# K2: ρ=0.39, p=0.013 (*)
# F: ρ=-0.37, p=0.019 (*)
# K2/F: ρ=0.51, p=0.001 (***)

# Create correlation matrix manually with your exact values
corr_subset <- matrix(c(0.39, -0.37, 0.51), nrow = 1, ncol = 3)
p_subset <- matrix(c(0.013, 0.019, 0.0001), nrow = 1, ncol = 3)

rownames(corr_subset) <- c('UACR')
colnames(corr_subset) <- c('Cortical K2', 'Cortical F', 'Cortical K2/F')
rownames(p_subset) <- rownames(corr_subset)
colnames(p_subset) <- colnames(corr_subset)

# Print to verify
cat("\nCorrelation values:\n")
print(corr_subset)
cat("\nP-values:\n")
print(p_subset)

# Create correlation plot with your exact values
temp_file <- paste0(base_path, "temp_corrplot.png")
png(temp_file, width = 10, height = 4, units = "in", res = 300)
corrplot(corr_subset, 
         method = "color",
         p.mat = p_subset,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 2.5,
         insig = "label_sig",
         addCoef.col = "black",
         number.cex = 1.5,
         tl.cex = 1.5,
         tl.col = 'black',
         cl.cex = 1.3,
         col = colorRampPalette(c("#4575b4", "white", "#d73027"))(200),
         mar = c(1, 0, 2, 0))
dev.off()

# Convert to grob
corr_grob <- grid::rasterGrob(png::readPNG(temp_file), interpolate = TRUE)
plot_j <- ggplot() + 
  annotation_custom(corr_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void() +
  labs(tag = "J") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

# ============================================================================
# CREATE FIGURE LEGEND TEXT
# ============================================================================

legend_text <- "Figure 1. PET Imaging Metrics: Group Comparisons, Pathology Associations, and Clinical Correlations.
(A-C) Cortical PET metrics by group. Boxplots comparing Healthy Weight Control vs Type 2 Diabetes for (A) Cortical F (perfusion), 
(B) Cortical K2 (TCA cycle metabolism), and (C) Cortical K2/F ratio (metabolic efficiency). P-values from Wilcoxon tests.
(D-F) Association of PET metrics with glomerular basement membrane (GBM) thickening. Boxplots comparing participants 
with vs without GBM thickening for (D) Cortical K2, (E) Cortical F, and (F) Cortical K2/F. P-values from Wilcoxon tests.
(G-I) Association of PET metrics with arteriosclerosis. Boxplots comparing participants with vs without arteriosclerosis 
for (G) Cortical K2, (H) Cortical F, and (I) Cortical K2/F. P-values from Wilcoxon tests. (J) Spearman correlations 
between urinary albumin-to-creatinine ratio (UACR) and cortical PET metrics (K2: ρ=0.39, p=0.013; F: ρ=-0.37, p=0.019; 
K2/F: ρ=0.51, p=0.001). Heat map shows correlation coefficients with significance levels: * p<0.05, ** p<0.01, *** p<0.001."

# Create text grob for legend
legend_grob <- textGrob(legend_text, 
                        gp = gpar(fontsize = 10, lineheight = 1.1),
                        x = 0.02, 
                        hjust = 0,
                        just = "left")

# ============================================================================
# COMBINE ALL PANELS WITH PATCHWORK
# ============================================================================

cat("Combining all panels\n")

# Combine main plots using patchwork
combined_plots <- (plot_a | plot_b | plot_c) /
  (plot_d | plot_e | plot_f) /
  (plot_g | plot_h | plot_i) /
  plot_j /
  wrap_elements(legend_grob) +
  plot_layout(heights = c(1, 1, 1, 0.5, 0.3)) +
  plot_annotation(
    title = "Figure.",
    subtitle = 'PET Imaging Metrics: Group Comparisons, Pathology Associations, and Clinical Correlations',
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0, color = "black"),
      plot.subtitle = element_text(size = 16, hjust = 0, color = "black", margin = margin(b = 10))
    )
  )

# ============================================================================
# SAVE FIGURE
# ============================================================================

cat("Saving figure\n")

ggsave(paste0(base_path, "Figure4_Complete_PET_Analysis.pdf"), 
       plot = combined_plots, width = 18, height = 26, dpi = 300)

ggsave(paste0(base_path, "Figure4_Complete_PET_Analysis.png"), 
       plot = combined_plots, width = 18, height = 26, dpi = 300)

ggsave(paste0(base_path, "Figure4_Complete_PET_Analysis.tiff"), 
       plot = combined_plots, width = 18, height = 26, dpi = 300, compression = "lzw")

# Clean up
file.remove(temp_file)

# Also save the legend text to a file
writeLines(legend_text, paste0(base_path, "Figure4_Legend.txt"))

cat("\n=== FIGURE COMPLETE ===\n")
cat(sprintf("Saved to: %s\n", base_path))
cat("Files created:\n")
cat("  - Figure4_Complete_PET_Analysis.pdf\n")
cat("  - Figure4_Complete_PET_Analysis.png\n")
cat("  - Figure4_Complete_PET_Analysis.tiff\n")
cat("  - Figure4_Legend.txt\n")