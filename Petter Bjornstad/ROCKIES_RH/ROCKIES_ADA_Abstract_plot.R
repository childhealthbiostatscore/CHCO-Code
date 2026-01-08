# ============================================================================
# MULTI-PANEL FIGURE: PET METRICS BY GROUP, PATHOLOGY, AND UACR
# Figure 1 - Updated with K2/K1 metric and reordered variables
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

# Calculate PET averages - UPDATED to include K1 and K2/K1
PET_avg <- function(data){
  tmp_df <- data %>% dplyr::select(lc_k2, rc_k2, lm_k2, rm_k2,
                                   lc_f, rc_f, lm_f, rm_f,
                                   lc_k1, rc_k1, lm_k1, rm_k1)
  avg_c_k2 <- tmp_df %>%
    dplyr::select(lc_k2, rc_k2) %>% rowMeans(na.rm=T)
  avg_m_k2 <- tmp_df %>% 
    dplyr::select(lm_k2, rm_k2) %>% rowMeans(na.rm=T)
  avg_c_f <- tmp_df %>% 
    dplyr::select(lc_f, rc_f) %>% rowMeans(na.rm=T)
  avg_m_f <- tmp_df %>% 
    dplyr::select(lm_f, rm_f) %>% rowMeans(na.rm=T)
  avg_c_k1 <- tmp_df %>% 
    dplyr::select(lc_k1, rc_k1) %>% rowMeans(na.rm=T)
  avg_m_k1 <- tmp_df %>% 
    dplyr::select(lm_k1, rm_k1) %>% rowMeans(na.rm=T)
  avg_c_k2_f <- avg_c_k2 / avg_c_f
  avg_m_k2_f <- avg_m_k2 / avg_m_f
  avg_c_k2_k1 <- avg_c_k2 / avg_c_k1
  avg_m_k2_k1 <- avg_m_k2 / avg_m_k1
  
  results <- bind_cols(avg_c_k2, avg_m_k2, avg_c_f, avg_m_f, 
                       avg_c_k1, avg_m_k1,
                       avg_c_k2_f, avg_m_k2_f,
                       avg_c_k2_k1, avg_m_k2_k1) %>% as.data.frame()
  names(results) <- c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 
                      'avg_c_k1', 'avg_m_k1',
                      'avg_c_k2_f', 'avg_m_k2_f',
                      'avg_c_k2_k1', 'avg_m_k2_k1')
  return(results)
}

tmp_results <- PET_avg(dat)

# Remove any existing PET average columns before binding new ones
dat_results <- dat %>% 
  select(-any_of(c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 
                   'avg_c_k1', 'avg_m_k1',
                   'avg_c_k2_f', 'avg_m_k2_f',
                   'avg_c_k2_k1', 'avg_m_k2_k1'))) %>%
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

# Calculate K1 averages from the FULL harmonized data (before SGLT2i filtering)
# This ensures all participants in dat_clinical get K2/K1 values
k1_data <- dat %>%
  filter(!is.na(lc_k1) | !is.na(rc_k1)) %>%
  mutate(avg_c_k1 = rowMeans(select(., lc_k1, rc_k1), na.rm = TRUE)) %>%
  select(record_id, avg_c_k1) %>%
  distinct(record_id, .keep_all = TRUE)

# Merge K1 into dat_clinical and calculate K2/K1
dat_clinical <- dat_clinical %>%
  left_join(k1_data, by = "record_id") %>%
  mutate(avg_c_k2_k1 = avg_c_k2 / avg_c_k1)

# Check how many now have K2/K1
cat(sprintf("N with K2/K1 in dat_clinical: %d / %d\n", 
            sum(!is.na(dat_clinical$avg_c_k2_k1)), nrow(dat_clinical)))

cat("All data loaded successfully\n\n")

# ============================================================================
# PANEL A-D: PET METRICS BY GROUP (Healthy Weight Control vs T2D)
# Order: F, K2, K2/F, K2/K1
# ============================================================================

cat("Creating Panels A-D: PET by Group\n")

# Prepare data for group comparison
aim1_df <- dat2 %>% 
  dplyr::select(record_id, group, avg_c_f, avg_c_k2, avg_c_k2_f, avg_c_k2_k1)

aim1_long <- aim1_df %>%
  pivot_longer(cols = c(avg_c_f, avg_c_k2, avg_c_k2_f, avg_c_k2_k1),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = factor(metric, 
                         levels = c("avg_c_f", "avg_c_k2", "avg_c_k2_f", "avg_c_k2_k1"),
                         labels = c("Cortical F", "Cortical K2", "Cortical K2/F", "Cortical K2/K1")))

# Create group comparison plots
group_colors <- c("Healthy Weight Control" = "#c2dfe3", "Type 2 Diabetes" = "#fcb1a6")

# Panel A: Cortical F
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

# Panel B: Cortical K2
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

# Panel C: Cortical K2/F
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

# Panel D: Cortical K2/K1 (NEW)
plot_d <- aim1_long %>% filter(metric == "Cortical K2/K1") %>%
  ggplot(aes(x = group, y = value, fill = group)) +
  geom_boxplot(color = "black", alpha = 0.8, outlier.shape = NA, linewidth = 0.5) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6, size = 2, shape = 21, color = "black", stroke = 0.3) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5) +
  labs(x = NULL, y = "Cortical K2/K1", tag = "D") +
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
# PANEL E-H: GBM THICKNESS (Order: F, K2, K2/F, K2/K1)
# ============================================================================

cat("Creating Panels E-H: GBM Thickness\n")

combined_df_gbm <- dat_clinical %>% 
  dplyr::select(record_id, avg_c_k2, avg_c_f, avg_c_k2_f, avg_c_k2_k1, `GBM thickness`) %>% 
  filter(`GBM thickness` != '', !is.na(`GBM thickness`))

pathology_colors <- c("no" = "#00008B", "yes" = "#8B0000")

# Panel E: Cortical F
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

# Panel F: Cortical K2
plot_f <- ggplot(combined_df_gbm, aes(x = `GBM thickness`, y = avg_c_k2, fill = `GBM thickness`)) +
  geom_boxplot(color = "black", alpha = 0.8, linewidth = 0.5) +
  scale_fill_manual(values = pathology_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5, p.adjust.method = "none") +
  labs(x = "GBM Thickening", y = "Cortical K2", fill = "GBM Thickening", tag = "F") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none",
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Panel G: Cortical K2/F
plot_g <- ggplot(combined_df_gbm, aes(x = `GBM thickness`, y = avg_c_k2_f, fill = `GBM thickness`)) +
  geom_boxplot(color = "black", alpha = 0.8, linewidth = 0.5) +
  scale_fill_manual(values = pathology_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5, p.adjust.method = "none") +
  labs(x = "GBM Thickening", y = "Cortical K2/F", fill = "GBM Thickening", tag = "G") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none",
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Panel H: Cortical K2/K1 (NEW)
plot_h <- ggplot(combined_df_gbm, aes(x = `GBM thickness`, y = avg_c_k2_k1, fill = `GBM thickness`)) +
  geom_boxplot(color = "black", alpha = 0.8, linewidth = 0.5) +
  scale_fill_manual(values = pathology_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5, p.adjust.method = "none") +
  labs(x = "GBM Thickening", y = "Cortical K2/K1", fill = "GBM Thickening", tag = "H") +
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
# PANEL I-L: ARTERIOSCLEROSIS (Order: F, K2, K2/F, K2/K1)
# ============================================================================

cat("Creating Panels I-L: Arteriosclerosis\n")

combined_df_arterio <- dat_clinical %>% 
  dplyr::select(record_id, avg_c_k2, avg_c_f, avg_c_k2_f, avg_c_k2_k1, arteriosclerosis) %>% 
  filter(arteriosclerosis != '', !is.na(arteriosclerosis))

# Panel I: Cortical F
plot_i <- ggplot(combined_df_arterio, aes(x = arteriosclerosis, y = avg_c_f, fill = arteriosclerosis)) +
  geom_boxplot(color = "black", alpha = 0.8, linewidth = 0.5) +
  scale_fill_manual(values = pathology_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5, p.adjust.method = "none") +
  labs(x = "Arteriosclerosis", y = "Cortical F", fill = "Arteriosclerosis", tag = "I") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none",
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Panel J: Cortical K2
plot_j <- ggplot(combined_df_arterio, aes(x = arteriosclerosis, y = avg_c_k2, fill = arteriosclerosis)) +
  geom_boxplot(color = "black", alpha = 0.8, linewidth = 0.5) +
  scale_fill_manual(values = pathology_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5, p.adjust.method = "none") +
  labs(x = "Arteriosclerosis", y = "Cortical K2", fill = "Arteriosclerosis", tag = "J") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none",
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Panel K: Cortical K2/F
plot_k <- ggplot(combined_df_arterio, aes(x = arteriosclerosis, y = avg_c_k2_f, fill = arteriosclerosis)) +
  geom_boxplot(color = "black", alpha = 0.8, linewidth = 0.5) +
  scale_fill_manual(values = pathology_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5, p.adjust.method = "none") +
  labs(x = "Arteriosclerosis", y = "Cortical K2/F", fill = "Arteriosclerosis", tag = "K") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none",
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Panel L: Cortical K2/K1 (NEW)
plot_l <- ggplot(combined_df_arterio, aes(x = arteriosclerosis, y = avg_c_k2_k1, fill = arteriosclerosis)) +
  geom_boxplot(color = "black", alpha = 0.8, linewidth = 0.5) +
  scale_fill_manual(values = pathology_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5, p.adjust.method = "none") +
  labs(x = "Arteriosclerosis", y = "Cortical K2/K1", fill = "Arteriosclerosis", tag = "L") +
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
# PANEL M: UACR CORRELATIONS (Order: F, K2, K2/F, K2/K1)
# Using ggplot heatmap instead of corrplot for better patchwork integration
# ============================================================================

cat("Creating Panel M: UACR Correlations\n")

# Calculate K2/K1 correlation with UACR from your data
uacr_k2k1_corr <- cor.test(dat2$acr_u, dat2$avg_c_k2_k1, method = "spearman", use = "complete.obs")

# Create data frame for ggplot heatmap
# Order: F, K2, K2/F, K2/K1
corr_df <- data.frame(
  metric = factor(c('Cortical F', 'Cortical K2', 'Cortical K2/F', 'Cortical K2/K1'),
                  levels = c('Cortical F', 'Cortical K2', 'Cortical K2/F', 'Cortical K2/K1')),
  rho = c(-0.37, 0.39, 0.51, round(uacr_k2k1_corr$estimate, 2)),
  pval = c(0.019, 0.013, 0.001, uacr_k2k1_corr$p.value)
)

# Add significance stars
corr_df <- corr_df %>%
  mutate(sig = case_when(
    pval < 0.001 ~ "***",
    pval < 0.01 ~ "**",
    pval < 0.05 ~ "*",
    TRUE ~ ""
  ),
  label = paste0(rho, sig))

# Print to verify
cat("\nCorrelation values:\n")
print(corr_df)

# Create ggplot heatmap
plot_m <- ggplot(corr_df, aes(x = metric, y = "UACR", fill = rho)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 6, fontface = "bold") +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", 
                       midpoint = 0, limits = c(-1, 1),
                       name = "Spearman ρ") +
  labs(x = NULL, y = NULL, tag = "M") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.tag = element_text(size = 16, face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# ============================================================================
# COMBINE ALL PANELS WITH PATCHWORK
# ============================================================================

cat("Combining all panels\n")

# Combine main plots using patchwork (now 4 columns per row)
combined_plots <- (plot_a | plot_b | plot_c | plot_d) /
  (plot_e | plot_f | plot_g | plot_h) /
  (plot_i | plot_j | plot_k | plot_l) /
  plot_m +
  plot_layout(heights = c(1, 1, 1, 0.5))

# ============================================================================
# SAVE FIGURE
# ============================================================================

cat("Saving figure\n")

ggsave(paste0(base_path, "Figure1_Complete_PET_Analysis.pdf"), 
       plot = combined_plots, width = 22, height = 26, dpi = 300)

ggsave(paste0(base_path, "Figure1_Complete_PET_Analysis.png"), 
       plot = combined_plots, width = 22, height = 26, dpi = 300)

ggsave(paste0(base_path, "Figure1_Complete_PET_Analysis.tiff"), 
       plot = combined_plots, width = 22, height = 26, dpi = 300, compression = "lzw")

# ============================================================================
# OUTPUT SUMMARY STATISTICS FOR ABSTRACT
# ============================================================================

cat("\n=== SUMMARY STATISTICS FOR ABSTRACT ===\n\n")

# Group comparisons (Healthy Weight Control vs T2D)
cat("--- GROUP COMPARISONS (Control vs T2D) ---\n\n")

# Calculate means by group
group_summary <- aim1_df %>%
  group_by(group) %>%
  summarise(
    n = n(),
    F_mean = mean(avg_c_f, na.rm = TRUE),
    F_sd = sd(avg_c_f, na.rm = TRUE),
    K2_mean = mean(avg_c_k2, na.rm = TRUE),
    K2_sd = sd(avg_c_k2, na.rm = TRUE),
    K2_F_mean = mean(avg_c_k2_f, na.rm = TRUE),
    K2_F_sd = sd(avg_c_k2_f, na.rm = TRUE),
    K2_K1_mean = mean(avg_c_k2_k1, na.rm = TRUE),
    K2_K1_sd = sd(avg_c_k2_k1, na.rm = TRUE)
  )

print(group_summary)

# Calculate percent differences (T2D vs Control)
control_vals <- group_summary %>% filter(group == "Healthy Weight Control")
t2d_vals <- group_summary %>% filter(group == "Type 2 Diabetes")

cat("\n--- PERCENT DIFFERENCES (T2D vs Control) ---\n")
cat(sprintf("F: %.1f%% difference\n", (t2d_vals$F_mean - control_vals$F_mean) / control_vals$F_mean * 100))
cat(sprintf("K2: %.1f%% difference\n", (t2d_vals$K2_mean - control_vals$K2_mean) / control_vals$K2_mean * 100))
cat(sprintf("K2/F: %.1f%% difference\n", (t2d_vals$K2_F_mean - control_vals$K2_F_mean) / control_vals$K2_F_mean * 100))
cat(sprintf("K2/K1: %.1f%% difference\n", (t2d_vals$K2_K1_mean - control_vals$K2_K1_mean) / control_vals$K2_K1_mean * 100))

# Wilcoxon tests for group comparisons
cat("\n--- P-VALUES (Wilcoxon tests) ---\n")
f_test <- wilcox.test(avg_c_f ~ group, data = aim1_df)
k2_test <- wilcox.test(avg_c_k2 ~ group, data = aim1_df)
k2f_test <- wilcox.test(avg_c_k2_f ~ group, data = aim1_df)
k2k1_test <- wilcox.test(avg_c_k2_k1 ~ group, data = aim1_df)

cat(sprintf("F: p = %.4g\n", f_test$p.value))
cat(sprintf("K2: p = %.4g\n", k2_test$p.value))
cat(sprintf("K2/F: p = %.4g\n", k2f_test$p.value))
cat(sprintf("K2/K1: p = %.4g\n", k2k1_test$p.value))

# GBM Thickening associations
cat("\n--- GBM THICKENING ASSOCIATIONS ---\n")
if(nrow(combined_df_gbm) > 0) {
  gbm_summary <- combined_df_gbm %>%
    group_by(`GBM thickness`) %>%
    summarise(
      n = n(),
      F_mean = mean(avg_c_f, na.rm = TRUE),
      K2_mean = mean(avg_c_k2, na.rm = TRUE),
      K2_F_mean = mean(avg_c_k2_f, na.rm = TRUE),
      K2_K1_mean = mean(avg_c_k2_k1, na.rm = TRUE)
    )
  print(gbm_summary)
  
  cat("\nP-values:\n")
  cat(sprintf("F: p = %.4g\n", wilcox.test(avg_c_f ~ `GBM thickness`, data = combined_df_gbm)$p.value))
  cat(sprintf("K2: p = %.4g\n", wilcox.test(avg_c_k2 ~ `GBM thickness`, data = combined_df_gbm)$p.value))
  cat(sprintf("K2/F: p = %.4g\n", wilcox.test(avg_c_k2_f ~ `GBM thickness`, data = combined_df_gbm)$p.value))
  cat(sprintf("K2/K1: p = %.4g\n", wilcox.test(avg_c_k2_k1 ~ `GBM thickness`, data = combined_df_gbm)$p.value))
}

# Arteriosclerosis associations
cat("\n--- ARTERIOSCLEROSIS ASSOCIATIONS ---\n")
if(nrow(combined_df_arterio) > 0) {
  arterio_summary <- combined_df_arterio %>%
    group_by(arteriosclerosis) %>%
    summarise(
      n = n(),
      F_mean = mean(avg_c_f, na.rm = TRUE),
      K2_mean = mean(avg_c_k2, na.rm = TRUE),
      K2_F_mean = mean(avg_c_k2_f, na.rm = TRUE),
      K2_K1_mean = mean(avg_c_k2_k1, na.rm = TRUE)
    )
  print(arterio_summary)
  
  cat("\nP-values:\n")
  cat(sprintf("F: p = %.4g\n", wilcox.test(avg_c_f ~ arteriosclerosis, data = combined_df_arterio)$p.value))
  cat(sprintf("K2: p = %.4g\n", wilcox.test(avg_c_k2 ~ arteriosclerosis, data = combined_df_arterio)$p.value))
  cat(sprintf("K2/F: p = %.4g\n", wilcox.test(avg_c_k2_f ~ arteriosclerosis, data = combined_df_arterio)$p.value))
  cat(sprintf("K2/K1: p = %.4g\n", wilcox.test(avg_c_k2_k1 ~ arteriosclerosis, data = combined_df_arterio)$p.value))
}

# UACR Correlations
cat("\n--- UACR CORRELATIONS (Spearman) ---\n")
uacr_f <- cor.test(dat2$acr_u, dat2$avg_c_f, method = "spearman", use = "complete.obs")
uacr_k2 <- cor.test(dat2$acr_u, dat2$avg_c_k2, method = "spearman", use = "complete.obs")
uacr_k2f <- cor.test(dat2$acr_u, dat2$avg_c_k2_f, method = "spearman", use = "complete.obs")
uacr_k2k1 <- cor.test(dat2$acr_u, dat2$avg_c_k2_k1, method = "spearman", use = "complete.obs")

cat(sprintf("F: rho = %.2f, p = %.4g\n", uacr_f$estimate, uacr_f$p.value))
cat(sprintf("K2: rho = %.2f, p = %.4g\n", uacr_k2$estimate, uacr_k2$p.value))
cat(sprintf("K2/F: rho = %.2f, p = %.4g\n", uacr_k2f$estimate, uacr_k2f$p.value))
cat(sprintf("K2/K1: rho = %.2f, p = %.4g\n", uacr_k2k1$estimate, uacr_k2k1$p.value))

cat("\n=== FIGURE COMPLETE ===\n")
cat(sprintf("Saved to: %s\n", base_path))
cat("Files created:\n")
cat("  - Figure1_Complete_PET_Analysis.pdf\n")
cat("  - Figure1_Complete_PET_Analysis.png\n")
cat("  - Figure1_Complete_PET_Analysis.tiff\n")