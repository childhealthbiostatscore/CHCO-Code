### ROCKIES Updates 9/15/25


library(scran)
library(future)
library(future.apply)
library(tidyverse)
library(colorspace)
library(patchwork)
library(ggdendro)
library(cowplot)
library(ggpubr)
library(rstatix)
library(arsenal)
library(Biobase)
library(msigdbr)
library(kableExtra)
library(knitr)
library(REDCapR)
library(data.table)
library(emmeans)
library(NMF)
library(pheatmap)
library(UpSetR)
library(enrichR)
library(WriteXLS)
library(SAVER)
library(readxl)
library(limma)
library(edgeR)
library(BiocGenerics)
library(GSEABase)
library(slingshot)
library(SingleCellExperiment)
library(MAST)
library(muscat)
library(scater)
library(Seurat)
library(jsonlite)
library(dplyr)
library(glmmTMB)
library(reshape2)
library(broom.mixed)
library(nebula)

















#Step 1: Hypermetabolism in T2D vs. HC (PET Analysis)



#harmonized_data <- read.csv("C:/Users/netio/Documents/Harmonized_data/harmonized_dataset.csv", na = '')




harmonized_data <- read.csv("C:/Users/netio/Documents/Harmonized_data/harmonized_dataset.csv", na = '')

#harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

#date_of_screen
#screen_date

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


PET_avg <- function(data){
  tmp_df <- data %>% dplyr::select(lc_k2_wo_cyst_vw, rc_k2_wo_cyst_vw, lm_k2_wo_cyst_vw, rm_k2_wo_cyst_vw,
                                   lc_f, rc_f, lm_f, rm_f)
  avg_c_k2 <- tmp_df %>%
    dplyr::select(lc_k2_wo_cyst_vw, rc_k2_wo_cyst_vw) %>% rowMeans(na.rm=T)
  
  avg_m_k2 <- tmp_df %>% 
    dplyr::select(lm_k2_wo_cyst_vw, rm_k2_wo_cyst_vw) %>% rowMeans(na.rm=T)
  
  avg_c_f <- tmp_df %>% 
    dplyr::select(lc_f, rc_f) %>% rowMeans(na.rm=T)
  
  avg_m_f <- tmp_df %>% 
    dplyr::select(lm_f, rm_f) %>% rowMeans(na.rm=T)
  
  avg_c_k2_f <- avg_c_k2 / avg_c_f
  
  avg_m_k2_f <- avg_m_k2/ avg_m_f
  
  results <- bind_cols(avg_c_k2, avg_m_k2, avg_c_f, avg_m_f, 
                       avg_c_k2_f, avg_m_k2_f) %>% as.data.frame()
  names(results) <- c('avg_c_k2_vw', 'avg_m_k2_vw', 'avg_c_f_vw', 'avg_m_f_vw', 
                      'avg_c_k2_f_vw', 'avg_m_k2_f_vw')
  
  return(results)
  
}


tmp_results_vw <- PET_avg(dat)


dat_results <- dat




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



dat_results <- dat_results %>% bind_cols(tmp_results, tmp_results_vw)


dat_results <- dat_results %>% filter(!is.na(avg_c_k2))

dat_results <- dat_results %>% filter(group %in% c('Lean Control', 'Type 2 Diabetes'))
dat_results$group2 <- NA

need_med_info <- dat_results %>% filter(is.na(group2))

dat2 <- dat_results

RH <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/RENALHEIR-SGLT2.csv')
names(RH) <- c('Subject', 'rep_instr', 'rep_inst', 'SGLT2')
RH2 <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/RenalHEIRitage-SGLT2Use.csv')
names(RH2) <- c('Subject', 'event', 'rep_instr', 'rep_inst', 'mrn', 'SGLT2', 'SGLT2_ever')
RH2 <- RH2 %>% filter(!is.na(mrn))
improve <- data.table::fread('C:/Users/netio/Downloads/IMPROVET2D-SGLT2i_DATA_LABELS_2025-08-25_0938.csv')
names(improve)[5] <- 'SGLT2'
names(improve)[1] <- 'record_id'

improve <- improve %>% filter(!is.na(SGLT2)) %>%
  filter(SGLT2 != '')

improve_small <- improve %>% filter(record_id %in% need_med_info$record_id)
RH_small <- RH %>% filter(Subject %in% need_med_info$record_id)
RH2_small <- RH2 %>% filter(mrn %in% need_med_info$mrn)

for(i in c(1:nrow(RH_small))){
  if(nrow(RH_small) == 0){
    next
  }
  if(RH_small$SGLT2[i] == 'No'){
    dat2$group2[which(dat2$record_id == RH_small$Subject[i])] <- 'T2D-No SGLTi2'
    dat2$epic_sglti2_1[which(dat2$record_id == RH_small$Subject[i])] <- 'No'
  }else if(RH_small$SGLT2[i] == 'Yes'){
    dat2$group2[which(dat2$record_id == RH_small$Subject[i])] <- 'T2D-SGLTi2'
    dat2$epic_sglti2_1[which(dat2$record_id == RH_small$Subject[i])] <- 'Yes'
  }else{
    next
  }
}

for(i in c(1:nrow(RH2_small))){
  if(nrow(RH2_small) == 0){
    next
  }
  if(RH2_small$SGLT2[i] == 'No'){
    dat2$group2[which(dat2$mrn == RH2_small$mrn[i])] <- 'T2D-No SGLTi2'
    dat2$epic_sglti2_1[which(dat2$mrn == RH2_small$mrn[i])] <- 'No'
  }else if(RH2_small$SGLT2[i] == 'Yes'){
    dat2$group2[which(dat2$mrn == RH2_small$mrn[i])] <- 'T2D-SGLTi2'
    dat2$epic_sglti2_1[which(dat2$mrn == RH2_small$mrn[i])] <- 'Yes'
  }else{
    next
  }
}

for(i in c(1:nrow(improve_small))){
  if(nrow(improve_small) == 0){
    next
  }
  if(improve_small$SGLT2[i] == 'No'){
    dat2$group2[which(dat2$record_id == improve_small$record_id[i])] <- 'T2D-No SGLTi2'
    dat2$epic_sglti2_1[which(dat2$record_id == improve_small$record_id[i])] <- 'No'
  }else if(improve_small$SGLT2[i] == 'Yes'){
    dat2$group2[which(dat2$record_id == improve_small$record_id[i])] <- 'T2D-SGLTi2'
    dat2$epic_sglti2_1[which(dat2$record_id == improve_small$record_id[i])] <- 'Yes'
  }else{
    next
  }
}


dat2$epic_sglti2_1[which(dat2$group == 'Lean Control')] <- 'No'

dat2 <- dat2 %>% filter(epic_sglti2_1 != 'Yes')


table1::table1(~age + sex + bmi +hba1c +  study + epic_sglti2_1 + avg_c_k2 + avg_c_k2_f | group, 
               data = dat2)



aim1_df <- dat2 %>% 
  dplyr::select(record_id, group, avg_c_f, avg_c_k2, avg_c_k2_f)

library(tidyverse)
library(ggpubr)
library(rstatix)

# Reshape data from wide to long format
aim1_long <- aim1_df %>%
  pivot_longer(cols = c(avg_c_f, avg_c_k2, avg_c_k2_f),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = factor(metric, 
                         levels = c("avg_c_f", "avg_c_k2", "avg_c_k2_f"),
                         labels = c("Cortical F", "Cortical K2", "Cortical K2/F")))

# Create base plot with boxplots and points
p <- ggplot(aim1_long, aes(x = metric, y = value)) +
  geom_boxplot(aes(fill = group), 
               color = "black", 
               alpha = 0.8, 
               outlier.shape = NA,
               position = position_dodge(width = 0.8)) +
  geom_point(aes(fill = group), 
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, 
             size = 1.5, 
             shape = 21, 
             color = "black") +
  scale_fill_manual(values = c("#c2dfe3", "#fcb1a6")) +
  labs(title = "Voxel-Based PET Imaging Metrics without Cysts by Group",
       x = "PET Metrics",
       y = "Value",
       fill = "Group") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# Calculate statistics for annotations
stat_test <- aim1_long %>%
  group_by(metric) %>%
  pairwise_wilcox_test(value ~ group, p.adjust.method = "none") %>%
  add_significance() %>%
  add_xy_position(x = "metric", dodge = 0.8)

# Adjust statistical test positions for the broken axis
stat_test_adjusted <- stat_test %>%
  mutate(y.position = case_when(
    metric == "Cortical F" ~ y.position,  # Keep original position for high values
    metric == "Cortical K2" ~ pmax(y.position, 0.20),       # Ensure minimum height in lower range
    metric == "Cortical K2/F" ~ pmax(0.20),     # Ensure minimum height in lower range
    TRUE ~ y.position
  ))

stat_test_adjusted$y.position <- c(3, 0.28, 0.30) 

# Create plot with statistical annotations
p_with_stats <- p + 
  stat_pvalue_manual(stat_test_adjusted, 
                     label = "p.adj.signif",
                     tip.length = 0.01,
                     step.increase = 0.005,
                     hide.ns = TRUE)

# Display the regular plot
print(p_with_stats)

# Create segmented y-axis plot
library(ggbreak)

p_broken <- p_with_stats + 
  scale_y_break(c(0.35, 0.75), scales = 2) +
  theme(axis.text.y = element_text(size = 10))

print(p_broken)

pdf('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.16.25/Aim1_VoxelPET.pdf')
print(p_broken)
dev.off()


png('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.16.25/Aim1_VoxelPET.png')
print(p_broken)
dev.off()


#PET global

aim1_df <- dat2 %>% 
  dplyr::select(record_id, group, avg_c_f = avg_c_f_vw, avg_c_k2 = avg_c_k2_vw, avg_c_k2_f = avg_c_k2_f_vw)

library(tidyverse)
library(ggpubr)
library(rstatix)

# Reshape data from wide to long format
aim1_long <- aim1_df %>%
  pivot_longer(cols = c(avg_c_f, avg_c_k2, avg_c_k2_f),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = factor(metric, 
                         levels = c("avg_c_f", "avg_c_k2", "avg_c_k2_f"),
                         labels = c("Cortical F", "Cortical K2", "Cortical K2/F")))

# Create base plot with boxplots and points
p <- ggplot(aim1_long, aes(x = metric, y = value)) +
  geom_boxplot(aes(fill = group), 
               color = "black", 
               alpha = 0.8, 
               outlier.shape = NA,
               position = position_dodge(width = 0.8)) +
  geom_point(aes(fill = group), 
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, 
             size = 1.5, 
             shape = 21, 
             color = "black") +
  scale_fill_manual(values = c("#c2dfe3",  "#fcb1a6")) +
  labs(title = "Global PET Imaging Metrics by Group",
       x = "PET Metrics",
       y = "Value",
       fill = "Group") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# Calculate statistics for annotations
stat_test <- aim1_long %>%
  group_by(metric) %>%
  pairwise_wilcox_test(value ~ group, p.adjust.method = "none") %>%
  add_significance() %>%
  add_xy_position(x = "metric", dodge = 0.8)

# Adjust statistical test positions for the broken axis
stat_test_adjusted <- stat_test %>%
  mutate(y.position = case_when(
    metric == "Cortical F" ~ y.position,  # Keep original position for high values
    metric == "Cortical K2" ~ 0.20,       # Ensure minimum height in lower range
    metric == "Cortical K2/F" ~ 0.20,     # Ensure minimum height in lower range
    TRUE ~ y.position
  ))


stat_test_adjusted$y.position <- c(2.9, 0.25, 0.28)

# Create plot with statistical annotations
p_with_stats <- p + 
  stat_pvalue_manual(stat_test_adjusted, 
                     label = "p.adj.signif",
                     tip.length = 0.01,
#                    step.increase = 0.02,
                     hide.ns = TRUE)

# Display the regular plot
print(p_with_stats)

# Create segmented y-axis plot
library(ggbreak)

p_broken <- p_with_stats + 
  scale_y_break(c(0.3, 0.70), scales = 2) +
  theme(axis.text.y = element_text(size = 10))

print(p_broken)

pdf('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.16.25/Aim1_GlobalPET.pdf')
print(p_broken)
dev.off()


png('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.16.25/Aim1_GlobalPET.png')
print(p_broken)
dev.off()


















#Step 2: uACR associations with PET metabolism (GBM thickness, atherosclerosis)

remove(list=ls())






harmonized_data <- read.csv("C:/Users/netio/Documents/Harmonized_data/harmonized_dataset.csv", na = '')

#harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

#date_of_screen
#screen_date

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


PET_avg <- function(data){
  tmp_df <- data %>% dplyr::select(lc_k2_wo_cyst_vw, rc_k2_wo_cyst_vw, lm_k2_wo_cyst_vw, rm_k2_wo_cyst_vw,
                                   lc_f, rc_f, lm_f, rm_f)
  avg_c_k2 <- tmp_df %>%
    dplyr::select(lc_k2_wo_cyst_vw, rc_k2_wo_cyst_vw) %>% rowMeans(na.rm=T)
  
  avg_m_k2 <- tmp_df %>% 
    dplyr::select(lm_k2_wo_cyst_vw, rm_k2_wo_cyst_vw) %>% rowMeans(na.rm=T)
  
  avg_c_f <- tmp_df %>% 
    dplyr::select(lc_f, rc_f) %>% rowMeans(na.rm=T)
  
  avg_m_f <- tmp_df %>% 
    dplyr::select(lm_f, rm_f) %>% rowMeans(na.rm=T)
  
  avg_c_k2_f <- avg_c_k2 / avg_c_f
  
  avg_m_k2_f <- avg_m_k2/ avg_m_f
  
  results <- bind_cols(avg_c_k2, avg_m_k2, avg_c_f, avg_m_f, 
                       avg_c_k2_f, avg_m_k2_f) %>% as.data.frame()
  names(results) <- c('avg_c_k2_vw', 'avg_m_k2_vw', 'avg_c_f_vw', 'avg_m_f_vw', 
                      'avg_c_k2_f_vw', 'avg_m_k2_f_vw')
  
  return(results)
  
}


tmp_results_vw <- PET_avg(dat)


dat_results <- dat




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



dat_results <- dat_results %>% bind_cols(tmp_results, tmp_results_vw)


dat_results <- dat_results %>% filter(!is.na(avg_c_k2))

dat_results <- dat_results %>% filter(group %in% c('Lean Control', 'Type 2 Diabetes'))
dat_results$group2 <- NA

need_med_info <- dat_results %>% filter(is.na(group2))

dat2 <- dat_results

RH <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/RENALHEIR-SGLT2.csv')
names(RH) <- c('Subject', 'rep_instr', 'rep_inst', 'SGLT2')
RH2 <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/RenalHEIRitage-SGLT2Use.csv')
names(RH2) <- c('Subject', 'event', 'rep_instr', 'rep_inst', 'mrn', 'SGLT2', 'SGLT2_ever')
RH2 <- RH2 %>% filter(!is.na(mrn))
improve <- data.table::fread('C:/Users/netio/Downloads/IMPROVET2D-SGLT2i_DATA_LABELS_2025-08-25_0938.csv')
names(improve)[5] <- 'SGLT2'
names(improve)[1] <- 'record_id'

improve <- improve %>% filter(!is.na(SGLT2)) %>%
filter(SGLT2 != '')

improve_small <- improve %>% filter(record_id %in% need_med_info$record_id)
RH_small <- RH %>% filter(Subject %in% need_med_info$record_id)
RH2_small <- RH2 %>% filter(mrn %in% need_med_info$mrn)

for(i in c(1:nrow(RH_small))){
  if(nrow(RH_small) == 0){
    next
  }
if(RH_small$SGLT2[i] == 'No'){
dat2$group2[which(dat2$record_id == RH_small$Subject[i])] <- 'T2D-No SGLTi2'
dat2$epic_sglti2_1[which(dat2$record_id == RH_small$Subject[i])] <- 'No'
}else if(RH_small$SGLT2[i] == 'Yes'){
dat2$group2[which(dat2$record_id == RH_small$Subject[i])] <- 'T2D-SGLTi2'
dat2$epic_sglti2_1[which(dat2$record_id == RH_small$Subject[i])] <- 'Yes'
}else{
next
}
}

for(i in c(1:nrow(RH2_small))){
  if(nrow(RH2_small) == 0){
    next
  }
if(RH2_small$SGLT2[i] == 'No'){
dat2$group2[which(dat2$mrn == RH2_small$mrn[i])] <- 'T2D-No SGLTi2'
dat2$epic_sglti2_1[which(dat2$mrn == RH2_small$mrn[i])] <- 'No'
}else if(RH2_small$SGLT2[i] == 'Yes'){
dat2$group2[which(dat2$mrn == RH2_small$mrn[i])] <- 'T2D-SGLTi2'
dat2$epic_sglti2_1[which(dat2$mrn == RH2_small$mrn[i])] <- 'Yes'
}else{
next
}
}

for(i in c(1:nrow(improve_small))){
  if(nrow(improve_small) == 0){
    next
  }
if(improve_small$SGLT2[i] == 'No'){
dat2$group2[which(dat2$record_id == improve_small$record_id[i])] <- 'T2D-No SGLTi2'
dat2$epic_sglti2_1[which(dat2$record_id == improve_small$record_id[i])] <- 'No'
}else if(improve_small$SGLT2[i] == 'Yes'){
dat2$group2[which(dat2$record_id == improve_small$record_id[i])] <- 'T2D-SGLTi2'
dat2$epic_sglti2_1[which(dat2$record_id == improve_small$record_id[i])] <- 'Yes'
}else{
next
}
}


dat2$epic_sglti2_1[which(dat2$group == 'Lean Control')] <- 'No'

dat2 <- dat2 %>% filter(epic_sglti2_1 != 'Yes')

table1::table1(~age + sex + bmi +hba1c +  study + epic_sglti2_1 + acr_u + gbm_thick_artmean + gbm_thick_harmmean | group, 
               data = dat2)


aim2_df <- dat2%>% 
  dplyr::select(record_id, group, 
                acr_u, gbm_thick_artmean, gbm_thick_harmmean,
                starts_with('avg_c_'))

library(corrplot)

# Function to calculate correlation p-values with error handling
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Check if both columns have enough finite observations
      valid_pairs <- complete.cases(mat[, i], mat[, j])
      if (sum(valid_pairs) >= 3) {  # Need at least 3 pairs for correlation test
        tryCatch({
          tmp <- cor.test(mat[, i], mat[, j], ...)
          p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }, error = function(e) {
          p.mat[i, j] <<- p.mat[j, i] <<- NA
        })
      } else {
        p.mat[i, j] <- p.mat[j, i] <- NA
      }
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Calculate full correlation matrix and p-values (including NA columns)
full_corr <- cor(aim2_df[, 3:ncol(aim2_df)], use = 'pairwise.complete.obs')
full_p_mat <- cor.mtest(aim2_df[, 3:ncol(aim2_df)])

# Subset to your desired rows and columns (this will include NAs where data is missing)
aim2_corr_df <- full_corr[c(1:3), c(10, 9, 11, 13, 12, 14)]
aim2_p_mat <- full_p_mat[c(1:3), c(10, 9, 11, 13, 12, 14)]

# Apply your labels
colnames(aim2_corr_df) <- c('Cortical F', 'Cortical K2', 'Cortical K2/F', 
                            'Cortical F (voxel)', 'Cortical K2 (voxel)', 'Cortical K2/F (voxel)')
rownames(aim2_corr_df) <- c('Urine Albumin-Creatinine Ratio', 'GBM Thickness (Arithmetic)', 'GBM Thickness (Harmonic)')

# Apply same labels to p-value matrix
colnames(aim2_p_mat) <- colnames(aim2_corr_df)
rownames(aim2_p_mat) <- rownames(aim2_corr_df)

# Create plot with "?" for missing data


pdf('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.16.25/Aim2_PETComparisons.pdf', width = 10, height = 14)
corrplot(aim2_corr_df, 
         method = "color", 
         tl.col = "black",
         tl.cex = 0.8,
         p.mat = aim2_p_mat,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 1.5,
         insig = "label_sig",
         pch.col = "black",
         na.label = "?",        # Show "?" for NA correlations
         na.label.col = "black") # Color of the "?" symbols        
dev.off()

png('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.16.25/Aim2_PETComparisons.png', width = 1000, height = 1400)
corrplot(aim2_corr_df, 
         method = "color", 
         tl.col = "black",
         tl.cex = 0.8,
         p.mat = aim2_p_mat,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 1.5,
         insig = "label_sig",
         pch.col = "black",
         na.label = "?",        # Show "?" for NA correlations
         na.label.col = "black") # Color of the "?" symbols        
dev.off()





#Step 3: Insulin clamp correlations with PET metabolism

remove(list=ls())




harmonized_data <- read.csv("C:/Users/netio/Documents/Harmonized_data/harmonized_dataset.csv", na = '')

#harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

#date_of_screen
#screen_date

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


PET_avg <- function(data){
  tmp_df <- data %>% dplyr::select(lc_k2_wo_cyst_vw, rc_k2_wo_cyst_vw, lm_k2_wo_cyst_vw, rm_k2_wo_cyst_vw,
                                   lc_f, rc_f, lm_f, rm_f)
  avg_c_k2 <- tmp_df %>%
    dplyr::select(lc_k2_wo_cyst_vw, rc_k2_wo_cyst_vw) %>% rowMeans(na.rm=T)
  
  avg_m_k2 <- tmp_df %>% 
    dplyr::select(lm_k2_wo_cyst_vw, rm_k2_wo_cyst_vw) %>% rowMeans(na.rm=T)
  
  avg_c_f <- tmp_df %>% 
    dplyr::select(lc_f, rc_f) %>% rowMeans(na.rm=T)
  
  avg_m_f <- tmp_df %>% 
    dplyr::select(lm_f, rm_f) %>% rowMeans(na.rm=T)
  
  avg_c_k2_f <- avg_c_k2 / avg_c_f
  
  avg_m_k2_f <- avg_m_k2/ avg_m_f
  
  results <- bind_cols(avg_c_k2, avg_m_k2, avg_c_f, avg_m_f, 
                       avg_c_k2_f, avg_m_k2_f) %>% as.data.frame()
  names(results) <- c('avg_c_k2_vw', 'avg_m_k2_vw', 'avg_c_f_vw', 'avg_m_f_vw', 
                      'avg_c_k2_f_vw', 'avg_m_k2_f_vw')
  
  return(results)
  
}


tmp_results_vw <- PET_avg(dat)


dat_results <- dat




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



dat_results <- dat_results %>% bind_cols(tmp_results, tmp_results_vw)


#dat_results <- dat_results %>% filter(!is.na(avg_c_k2))

dat_results <- dat_results %>% filter(group %in% c('Lean Control', 'Obese Control', 'Type 2 Diabetes'))


aim3_df <- dat_results %>% 
  dplyr::select(record_id, group, 
                acprg, airg, #raw_m, 
                p1_gc_leanm, p1_gc_m, p1_raw_leanm, p1_raw_m,
                p2_gc_leanm, p2_gc_m, p2_raw_leanm, p2_raw_m,
                starts_with('avg_c_'))

aim3_df <- aim3_df %>% mutate(p1_DI_leanm = p1_gc_leanm * airg, 
                              p1_DI = p1_gc_m * airg, 
                              p1_DI_raw_leanm = p1_raw_leanm * airg, 
                              p1_DI_raw = p1_raw_m * airg, 
                              p2_DI_leanm = p2_gc_leanm * airg, 
                              p2_DI = p2_gc_m * airg, 
                              p2_DI_raw_leanm = p2_raw_leanm * airg, 
                              p2_DI_raw = p2_raw_m * airg)


library(corrplot)

# Function to calculate correlation p-values with error handling
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Check if both columns have enough finite observations
      valid_pairs <- complete.cases(mat[, i], mat[, j])
      if (sum(valid_pairs) >= 3) {  # Need at least 3 pairs for correlation test
        tryCatch({
          tmp <- cor.test(mat[, i], mat[, j], ...)
          p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }, error = function(e) {
          p.mat[i, j] <<- p.mat[j, i] <<- NA
        })
      } else {
        p.mat[i, j] <- p.mat[j, i] <- NA
      }
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Calculate full correlation matrix and p-values (including NA columns)
full_corr <- cor(aim3_df[, 3:ncol(aim3_df)], use = 'pairwise.complete.obs')
full_p_mat <- cor.mtest(aim3_df[, 3:ncol(aim3_df)])

# Subset to your desired rows and columns (this will include NAs where data is missing)
aim3_corr_df <- full_corr[c(1:10), c(17, 16, 18, 20, 19, 21)]
aim3_p_mat <- full_p_mat[c(1:10), c(17, 16, 18, 20, 19, 21)]

# Apply your labels
colnames(aim3_corr_df) <- c('Cortical F', 'Cortical K2', 'Cortical K2/F', 
                            'Cortical F (voxel)', 'Cortical K2 (voxel)', 'Cortical K2/F (voxel)')
rownames(aim3_corr_df) <- c('Acute C-Peptide Response to Glucose', 'Acute Insulin Response to Glucose', 
                            'Phase 1 GC Lean M-Value', 'Phase 1 GC M-Value', 'Phase 1 Raw Lean M-Value', 'Phase 1 Raw M-Value', 
                            'Phase 2 GC Lean M-Value', 'Phase 2 GC M-Value', 'Phase 2 Raw Lean M-Value', 'Phase 2 Raw M-Value')

# Apply same labels to p-value matrix
colnames(aim3_p_mat) <- colnames(aim3_corr_df)
rownames(aim3_p_mat) <- rownames(aim3_corr_df)

# Create plot with "?" for missing data


pdf('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.16.25/Aim3_PETComparisons.pdf', width = 10, height = 14)
corrplot(aim3_corr_df, 
         method = "color", 
         tl.col = "black",
         tl.cex = 0.8,
         p.mat = aim3_p_mat,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 1.5,
         insig = "label_sig",
         pch.col = "black",
         na.label = "?",        # Show "?" for NA correlations
         na.label.col = "black") # Color of the "?" symbols        
dev.off()

png('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.16.25/Aim3_PETComparisons.png', width = 1000, height = 1400)
corrplot(aim3_corr_df, 
         method = "color", 
         tl.col = "black",
         tl.cex = 0.8,
         p.mat = aim3_p_mat,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 1.5,
         insig = "label_sig",
         pch.col = "black",
         na.label = "?",        # Show "?" for NA correlations
         na.label.col = "black") # Color of the "?" symbols        
dev.off()













#Step 4: TCA & OxPhos transcript and pathways in PT cells 

#Step 5: Gene scores, pathways scores associated with PET variables

#Step 6: ROCKIES correlations with OGIS, HOMA, tubular sodium with K2 and K2/F

#Step 7: SGLT2 inhibitor use with K2 and K2/F

remove(list = ls())

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


PET_avg <- function(data){
  tmp_df <- data %>% dplyr::select(lc_k2_wo_cyst_vw, rc_k2_wo_cyst_vw, lm_k2_wo_cyst_vw, rm_k2_wo_cyst_vw,
                                   lc_f, rc_f, lm_f, rm_f)
  avg_c_k2 <- tmp_df %>%
    dplyr::select(lc_k2_wo_cyst_vw, rc_k2_wo_cyst_vw) %>% rowMeans(na.rm=T)
  
  avg_m_k2 <- tmp_df %>% 
    dplyr::select(lm_k2_wo_cyst_vw, rm_k2_wo_cyst_vw) %>% rowMeans(na.rm=T)
  
  avg_c_f <- tmp_df %>% 
    dplyr::select(lc_f, rc_f) %>% rowMeans(na.rm=T)
  
  avg_m_f <- tmp_df %>% 
    dplyr::select(lm_f, rm_f) %>% rowMeans(na.rm=T)
  
  avg_c_k2_f <- avg_c_k2 / avg_c_f
  
  avg_m_k2_f <- avg_m_k2/ avg_m_f
  
  results <- bind_cols(avg_c_k2, avg_m_k2, avg_c_f, avg_m_f, 
                       avg_c_k2_f, avg_m_k2_f) %>% as.data.frame()
  names(results) <- c('avg_c_k2_vw', 'avg_m_k2_vw', 'avg_c_f_vw', 'avg_m_f_vw', 
                      'avg_c_k2_f_vw', 'avg_m_k2_f_vw')
  
  return(results)
  
}


tmp_results_vw <- PET_avg(dat)


dat_results <- dat




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



dat_results <- dat_results %>% bind_cols(tmp_results, tmp_results_vw)


dat_results <- dat_results %>% filter(!is.na(avg_c_k2))

dat_results <- dat_results %>% filter(group %in% c('Lean Control', 'Obese Control', 'Type 2 Diabetes'))





#Step 8: Modifiers of changes in K2 (OGIS, Tubular sodium, urine metabolites, fumarate/malate)

remove(list=ls())

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))

#fractional excretion of sodium is FeNA (%) = (Urine Sodium/Plasma Sodium)/(Urine Creatinine/Plasma Creatinine) x 100
#OGIS uses an excel calculator, need to do more research on this 

aim8_df <- dat %>% 
  dplyr::select(record_id, group, 
                malate, fumarate, homa_ir, eGFR_bedside_Schwartz,
                eGFR_CKD_epi,
                eGFR_fas_cr,
                eGFR_fas_cr_cysc,
                eGFR_Zap,
                eGFR_Schwartz,
                sodium_s, sodium_u, creat_u, creat_s
  )









