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




#harmonized_data <- read.csv("", na = '')

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

#date_of_screen
#screen_date

#dat <- harmonized_data %>% dplyr::select(-dob) %>% 
#  arrange(date_of_screen) %>% 
#  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
#                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
#                   .by = c(record_id, visit))



dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
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






# Fix data types before creating the table
library(gtsummary)
library(gt)
library(dplyr)

# Convert variables to proper data types
combined_df <- dat2 %>%
  mutate(
    # Ensure continuous variables are numeric
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    hba1c = as.numeric(hba1c),
    
    # Ensure categorical variables are factors or characters
    sex = as.factor(sex),
    race_ethnicity = as.factor(race_ethnicity),
    study = as.factor(study),
    group = as.factor(group),
    epic_sglti2_1 = as.factor(epic_sglti2_1)
  )



# Now create the table with proper data types
desc_table1_fixed <- combined_df %>%
  select(age, sex, race_ethnicity, bmi, hba1c, study, group, epic_sglti2_1) %>%
  tbl_summary(
    by = group,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      hba1c ~ "continuous",
      sex ~ "categorical",
      race_ethnicity ~ "categorical",
      study ~ "categorical",
      epic_sglti2_1 ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      sex ~ "Sex", 
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      study ~ "Study",
      epic_sglti2_1 ~ "SGLT2 Inhibitor Use"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "t.test"
    # Skip categorical p-values if they cause issues
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Group**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

# Save version with epic
desc_table1_fixed %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/demographics_aim1_with_epic_final.png", 
         vwidth = 1200, vheight = 800)










#table1::table1(~age + sex + bmi +hba1c +  study + epic_sglti2_1 + avg_c_k2 + avg_c_k2_f | group, 
#               data = dat2)


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

pdf('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/Aim1_VoxelPET.pdf')
print(p_broken)
dev.off()


png('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/Aim1_VoxelPET.png')
print(p_broken)
dev.off()


#PET global

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

pdf('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/Aim1_GlobalPET.pdf')
print(p_broken)
dev.off()


png('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/Aim1_GlobalPET.png')
print(p_broken)
dev.off()












##################################################################################################################################




#Step 4: TCA & OxPhos transcript and pathways in PT cells 

library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)



#barplots

lc_files <- list.files('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', pattern='csv')



t2d_files <- list.files('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2/', pattern = 'csv')


celltypes <- c('PT', 'PT-S1/S2', 'PT-S3', 'aPT')


overall_summary <- data.frame(celltype = celltypes)
overall_summary$TCA_total <- 27
overall_summary$OxPhos_total <- 10
overall_summary$TCA_flipped <- NA
overall_summary$TCA_sig <- NA
overall_summary$OxPhos_flipped <- NA
overall_summary$OxPhos_sig <- NA


for(i in c(1:length(celltypes))){
  
  celltype <- celltypes[i]
  if(celltype == 'PT'){
    text_size <- 12
    legend_position <- 'right'
    axis_angle <- 0
    hjust_val <- 0
  }else{
    text_size <- 15
    legend_position <- 'none'
    axis_angle <- 45
    hjust_val <- 1
  }
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  tmp_lc <- lc_files[str_which(lc_files, pattern = paste0('cycle_', celltype2, '_cells'))]
  tmp_t2d <- t2d_files[str_which(t2d_files, pattern = paste0('cycle_', celltype2, '_cells'))]
  
  #TCA
  tmp_lc_tca <- tmp_lc[str_which(tmp_lc, pattern = 'TCA')]
  tmp_t2d_tca <- tmp_t2d[str_which(tmp_t2d, pattern = 'TCA')]
  
  tca_lc <- data.table::fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', tmp_lc_tca))
  tca_t2d <- data.table::fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2/', tmp_t2d_tca))
  
  tca_lc <- tca_lc %>% dplyr::select(gene, logFC_lc = logFC_groupType_2_Diabetes, 
                                     pvalue_lc =  any_of(c("p_groupType_2_Diabetes", "pvalue")))
  
  tca_t2d <- tca_t2d %>% dplyr::select(gene, logFC_t2d = logFC_epic_sglti2_1Yes, 
                                       pvalue_t2d =  any_of(c("p_epic_sglti2_1Yes", "pvalue")))
  
  tca_full <- tca_lc %>% left_join(tca_t2d)
  
  tca_full_sig <- tca_full %>% 
    mutate(direction = ifelse(logFC_lc < 0, 'Down', 
                              ifelse(logFC_lc > 0, 'Up', NA)))
  
  # significant_traits <- tca_full %>% filter(pvalue_lc < 0.05 | pvalue_t2d < 0.05)
  #  opposite_1 <- tca_full %>% filter(logFC_lc > 0 & logFC_t2d < 0)
  #  opposite_2 <- tca_full %>% filter(logFC_lc < 0 & logFC_t2d > 0)
  
  #  tca_full_sig <- bind_rows(list(significant_traits, opposite_1, opposite_2)) %>% 
  #    filter(!duplicated(gene))
  
  #  overall_summary$TCA_flipped[i] <- tca_full_sig %>% filter((logFC_lc > 0 & logFC_t2d < 0) | (logFC_lc < 0 & logFC_t2d > 0)) %>% 
  #    nrow()
  #  overall_summary$TCA_sig[i] <- tca_full_sig %>% filter(pvalue_lc < 0.05 | pvalue_t2d < 0.05) %>% nrow()
  
  
  if(nrow(tca_full_sig) > 0){
    plot_df <- tca_full_sig %>% 
      gather(key = 'metric_condition', value = 'value', -gene, -direction) %>% 
      separate(metric_condition, into = c('metric', 'condition'), sep = '_') %>% 
      spread(key = metric, value = value) %>% mutate(Significance = ifelse(pvalue < 0.05, '*', ''))
    
    plot_df <- plot_df %>% filter(condition =='lc')
    
    tmp_plot <- ggplot(plot_df, aes(x=gene, y=logFC, fill=direction))+
      geom_bar(stat='identity', position='dodge')+
      geom_text(aes(label = Significance), 
                position = position_dodge(width = 0.9),
                vjust = -0.5, 
                size = 4)+
      labs(x='Gene', y='LogFC', title = paste0('TCA Gene Comparison in ', celltype2, ' Cells'))+
      scale_fill_manual(name = 'Comparison',
                        labels = c('Up' = 'Up-regulated', 'Down' = 'Down-regulated'),
                        values = c('Down' = '#d73027', 'Up' = '#4575b4'))+
      theme_classic()+theme(text = element_text(size = text_size), legend.position = legend_position, axis.text.x = element_text(angle = axis_angle, hjust = hjust_val))
    
    pdf(paste0('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/barplots/TCA_', celltype2, '_barplot.pdf'), 
        width = 12, height = 8)
    print(tmp_plot)
    dev.off()
    
    png(paste0('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/barplots/TCA_', celltype2, '_barplot.png'), 
        width = 1200, height = 800)
    print(tmp_plot)
    dev.off()
    
    
    
  }
  
  
  
  #OxPhos
  tmp_lc_oxphos <- tmp_lc[str_which(tmp_lc, pattern = 'PHOS_')]
  tmp_t2d_oxphos <- tmp_t2d[str_which(tmp_t2d, pattern = 'PHOS_')]
  
  oxphos_lc <- data.table::fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', tmp_lc_oxphos))
  oxphos_t2d <- data.table::fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2/', tmp_t2d_oxphos))
  
  oxphos_lc <- oxphos_lc %>% dplyr::select(gene, logFC_lc = logFC_groupType_2_Diabetes, 
                                           pvalue_lc = any_of(c("p_groupType_2_Diabetes", "pvalue")))
  
  oxphos_t2d <- oxphos_t2d %>% dplyr::select(gene, logFC_t2d = logFC_epic_sglti2_1Yes, 
                                             pvalue_t2d = any_of(c("p_epic_sglti2_1Yes", "pvalue")))
  
  oxphos_full <- oxphos_lc %>% left_join(oxphos_t2d)
  
  # significant_traits <- oxphos_full %>% filter(pvalue_lc < 0.05 | pvalue_t2d < 0.05)
  #  opposite_1 <- oxphos_full %>% filter(logFC_lc > 0 & logFC_t2d < 0)
  #  opposite_2 <- oxphos_full %>% filter(logFC_lc < 0 & logFC_t2d > 0)
  
  # oxphos_full_sig <- bind_rows(list(significant_traits, opposite_1, opposite_2)) %>% 
  #  filter(!duplicated(gene))
  
  
  #overall_summary$OxPhos_flipped[i] <- oxphos_full_sig %>% filter((logFC_lc > 0 & logFC_t2d < 0) | (logFC_lc < 0 & logFC_t2d > 0)) %>% 
  # nrow()
  #overall_summary$OxPhos_sig[i] <- oxphos_full_sig %>% filter(pvalue_lc < 0.05 | pvalue_t2d < 0.05) %>% nrow()
  
  oxphos_full_sig <- oxphos_full %>% 
    mutate(direction = ifelse(logFC_lc < 0, 'Down', 
                              ifelse(logFC_lc > 0, 'Up', NA)))
  
  if(nrow(oxphos_full_sig) > 0){
    plot_df <- oxphos_full_sig %>% 
      gather(key = 'metric_condition', value = 'value', -gene, -direction) %>% 
      separate(metric_condition, into = c('metric', 'condition'), sep = '_') %>% 
      spread(key = metric, value = value) %>% mutate(Significance = ifelse(pvalue < 0.05, '*', ''))
    
    plot_df <- plot_df %>% filter(condition =='lc')
    
    tmp_plot <- ggplot(plot_df, aes(x=gene, y=logFC, fill=direction))+
      geom_bar(stat='identity', position='dodge')+
      geom_text(aes(label = Significance), 
                position = position_dodge(width = 0.9),
                vjust = -0.5, 
                size = 4)+
      labs(x='Gene', y='LogFC', title = paste0('OxPhos Gene Comparison in ', celltype2, ' Cells'))+
      scale_fill_manual(name = 'Comparison',
                        labels = c('Up' = 'Up-regulated', 'Down' = 'Down-regulated'),
                        values = c('Down' = '#d73027', 'Up' = '#4575b4'))+
      theme_classic()+theme(text = element_text(size = text_size), legend.position = legend_position, axis.text.x = element_text(angle = axis_angle, hjust = hjust_val))
    
    pdf(paste0('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/barplots/oxphos_', celltype2, '_barplot.pdf'), 
        width = 12, height = 8)
    print(tmp_plot)
    dev.off()
    
    png(paste0('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/barplots/oxphos_', celltype2, '_barplot.png'), 
        width = 1200, height = 800)
    print(tmp_plot)
    dev.off()
    
    
    
    
  }
  
  
  
  
  print(paste0(celltype2, ' is done.'))
  
}


library(enrichplot)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)  # Human annotation, change for other species
library(DOSE)

library(ggplot2)
library(dplyr)

# Alternative libraries for IPA-style analysis
library(msigdbr)  # For MSigDB gene sets
library(fgsea)    # Fast GSEA implementation

library(dplyr)
library(stringr)

library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)
library(msigdbr) # For Hallmark gene sets
library(enrichplot)

  # ===============================
  # CREATE DOTPLOT FUNCTION
  # ===============================
  create_dotplot <- function(enrich_obj, plot_title, analysis_type = "ORA") {
    if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
      cat("No significant results for", plot_title, "\n")
      return(NULL)
    }
    
    df <- as.data.frame(enrich_obj)
    
    if (analysis_type == "ORA") {
      plot <- df %>%
        slice_head(n = 15) %>%
        mutate(
          Description = str_wrap(Description, width = 50),
          Description = str_remove(Description, "^HALLMARK_|^GOBP_|^REACTOME_"), # Clean pathway names
          Description = forcats::fct_reorder(Description, Count),
          log_padj = -log10(pvalue) # Using raw p-values instead of adjusted
        ) %>%
        ggplot(aes(x = Count, y = Description)) +
        geom_point(aes(size = Count, color = log_padj)) +
        scale_color_gradient(low = "lightblue", high = "darkred", 
                             name = "-log10(Adj. P)") +
        scale_size_continuous(name = "Gene Count", range = c(2, 10)) +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 8, color = "black"),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.title = element_text(color = "black"),
          plot.title = element_text(size = 12, hjust = 0.5, color = "black"),
          legend.position = "right",
          legend.text = element_text(color = "black"),
          legend.title = element_text(color = "black"),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_line(color = "grey90", size = 0.5),
          panel.grid.minor = element_line(color = "grey95", size = 0.3)
        ) +
        labs(
          title = plot_title,
          x = "Gene Count",
          y = "Pathway"
        )
    } else { # GSEA
      plot <- df %>%
        slice_head(n = 15) %>%
        mutate(
          Description = str_wrap(Description, width = 50),
          Description = str_remove(Description, "^HALLMARK_|^GOBP_|^REACTOME_"),
          Description = forcats::fct_reorder(Description, NES),
          log_padj = -log10(p.adjust),
          Direction = ifelse(NES > 0, "Upregulated", "Downregulated")
        ) %>%
        ggplot(aes(x = NES, y = Description)) +
        geom_point(aes(size = abs(NES), color = log_padj)) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.7, color = "gray50") +
        scale_color_gradient(low = "lightblue", high = "darkred", 
                             name = "-log10(P-value)") +
        scale_size_continuous(name = "Abs(NES)", range = c(2, 10)) +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 8, color = "black"),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.title = element_text(color = "black"),
          plot.title = element_text(size = 12, hjust = 0.5, color = "black"),
          legend.position = "right",
          legend.text = element_text(color = "black"),
          legend.title = element_text(color = "black"),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_line(color = "grey90", size = 0.5),
          panel.grid.minor = element_line(color = "grey95", size = 0.3)
        ) +
        labs(
          title = plot_title,
          x = "Normalized Enrichment Score (NES)",
          y = "Pathway"
        )
    }
    
    return(plot)
  }
  
  # ===============================
  # GENERATE ALL PLOTS
  # ===============================
 
#plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
#plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)

#file.copy(from=plots.png.paths, to="C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/NEBULA/")



########################################################## Combine PNG files 


library(cowplot)
library(magick)
library(ggplot2)

# Read images - including the main overview plot
img_main <- image_read("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/barplots/oxphos_PT_barplot.png")
img1 <- image_read("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/barplots/oxphos_PT_S1_S2_barplot.png")
img2 <- image_read("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/barplots/oxphos_PT_S3_barplot.png")
img3 <- image_read("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/barplots/oxphos_aPT_barplot.png")

# Convert to ggplot objects
p_main <- ggdraw() + draw_image(img_main)
p1 <- ggdraw() + draw_image(img1)
p2 <- ggdraw() + draw_image(img2)
p3 <- ggdraw() + draw_image(img3)

# Create bottom row with three subplots
bottom_row <- plot_grid(p1, p2, p3, 
                        ncol = 3, 
                        labels = c("B", "C", "D"), 
                        label_size = 14,
                        rel_widths = c(1, 1, 1))

# Combine main plot on top with bottom row
final_plot <- plot_grid(p_main, bottom_row,
                        ncol = 1,           # 1 column (stack vertically)
                        nrow = 2,           # 2 rows
                        labels = c("A", ""), # Label only the main plot
                        label_size = 16,
                        rel_heights = c(1, 1.2))  # Bottom row slightly taller

# Save at high resolution
ggsave("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/barplots/oxphos_PT_complete_layout.png",
       plot = final_plot,
       width = 15,          # Width in inches
       height = 12,         # Height in inches
       dpi = 300,           # High resolution
       units = "in")



# Read images - including the main overview plot
img_main <- image_read("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/barplots/TCA_PT_barplot.png")
img1 <- image_read("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/barplots/TCA_PT_S1_S2_barplot.png")
img2 <- image_read("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/barplots/TCA_PT_S3_barplot.png")
img3 <- image_read("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/barplots/TCA_aPT_barplot.png")

# Convert to ggplot objects
p_main <- ggdraw() + draw_image(img_main)
p1 <- ggdraw() + draw_image(img1)
p2 <- ggdraw() + draw_image(img2)
p3 <- ggdraw() + draw_image(img3)

# Create bottom row with three subplots
bottom_row <- plot_grid(p1, p2, p3, 
                        ncol = 3, 
                        labels = c("B", "C", "D"), 
                        label_size = 14,
                        rel_widths = c(1, 1, 1))

# Combine main plot on top with bottom row
final_plot <- plot_grid(p_main, bottom_row,
                        ncol = 1,           # 1 column (stack vertically)
                        nrow = 2,           # 2 rows
                        labels = c("A", ""), # Label only the main plot
                        label_size = 16,
                        rel_heights = c(1, 1.2))  # Bottom row slightly taller

# Save at high resolution
ggsave("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/barplots/TCA_PT_complete_layout.png",
       plot = final_plot,
       width = 15,          # Width in inches
       height = 12,         # Height in inches
       dpi = 300,           # High resolution
       units = "in")





















###########################################################################################################################



remove(list=ls())



#Step 5: Gene scores, pathways scores associated with PET variables

### Module Score Analysis
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
library(doParallel)



load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')

so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

#dat_groups <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_GroupAssignments.txt')
#dat_groups <- dat_groups %>% filter(group2 %in% c('Lean Control', 'T2D-No SGLTi2'))

#so_subset <- subset(so_subset, record_id == dat_groups$record_id)
test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))

#load('C:/Users/netio/Downloads/TCA_genes.txt')
#load('C:/Users/netio/Downloads/OxPhos_genes.txt')


dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/'





#Make sure exposure/independent/x variable or group variable is a factor variable
so_subset$group <- factor(so_subset$group)
#Make sure to set reference level
so_subset$group  <- relevel(so_subset$group ,ref="Lean_Control")

# Load required libraries
library(msigdbr)
library(Seurat)

# Your original GO terms (keep if you want both analyses)
insulin_sensitivity_go_terms <- c(
  "GO:0032868",  # response to insulin
  "GO:0008286",  # insulin receptor signaling pathway
  "GO:0046627",  # negative regulation of insulin receptor signaling pathway
  "GO:0046628",  # positive regulation of insulin receptor signaling pathway
  "GO:0005159"   # insulin-like growth factor receptor binding
)

# Kidney-relevant HALLMARK pathways
kidney_hallmark_terms <- c(
  # Metabolism-related
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_FATTY_ACID_METABOLISM", 
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_ADIPOGENESIS",
  "HALLMARK_PEROXISOME",
  
  # Immune/Inflammatory
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_ALLOGRAFT_REJECTION",
  
  # Cross-cutting (Metabolism-Immune)
  "HALLMARK_HYPOXIA",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
)

# Get HALLMARK gene sets from MSigDB
library(dplyr)

# Get HALLMARK gene sets from MSigDB
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")

# Extract genes properly using dplyr
kidney_hallmark_genes <- list()
for(pathway in kidney_hallmark_terms) {
  genes <- hallmark_sets %>%
    filter(gs_name == pathway) %>%
    pull(gene_symbol)
  
  kidney_hallmark_genes[[pathway]] <- unique(genes)
  
  cat("HALLMARK pathway:", pathway, "- Genes:", length(kidney_hallmark_genes[[pathway]]), "\n")
}



# Optional: Keep your GO analysis
library(org.Hs.eg.db)
library(AnnotationDbi)
insulin_sensitivity_genes <- list()
for(go_term in insulin_sensitivity_go_terms) {
  genes <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = go_term,
    columns = c("SYMBOL", "ENTREZID", "ENSEMBL"),
    keytype = "GOALL"
  )
  
  genes_clean <- genes[!is.na(genes$SYMBOL) & genes$SYMBOL != "", ]
  insulin_sensitivity_genes[[go_term]] <- unique(genes_clean$SYMBOL)
  
  cat("GO term:", go_term, "- Genes:", length(insulin_sensitivity_genes[[go_term]]), "\n")
}

# Add your existing TCA and OxPhos scores
so_subset <- AddModuleScore(object = so_subset, 
                            features = list(tca_genes), 
                            name = 'TCA_score')
so_subset <- AddModuleScore(object = so_subset, 
                            features = list(ox_phos_genes),
                            name = 'OxPhos_score')

# Create clean names for HALLMARK pathways
hallmark_names <- c(
  'HALLMARK_OXIDATIVE_PHOSPHORYLATION' = 'OxPhos_Hallmark',
  'HALLMARK_FATTY_ACID_METABOLISM' = 'Fatty_Acid_Metabolism',
  'HALLMARK_GLYCOLYSIS' = 'Glycolysis',
  'HALLMARK_MTORC1_SIGNALING' = 'mTORC1_Signaling',
  'HALLMARK_ADIPOGENESIS' = 'Adipogenesis',
  'HALLMARK_PEROXISOME' = 'Peroxisome',
  'HALLMARK_INFLAMMATORY_RESPONSE' = 'Inflammatory_Response',
  'HALLMARK_INTERFERON_GAMMA_RESPONSE' = 'IFN_Gamma_Response',
  'HALLMARK_INTERFERON_ALPHA_RESPONSE' = 'IFN_Alpha_Response',
  'HALLMARK_IL6_JAK_STAT3_SIGNALING' = 'IL6_JAK_STAT3',
  'HALLMARK_TNFA_SIGNALING_VIA_NFKB' = 'TNFa_NF_kB',
  'HALLMARK_COMPLEMENT' = 'Complement',
  'HALLMARK_ALLOGRAFT_REJECTION' = 'Allograft_Rejection',
  'HALLMARK_HYPOXIA' = 'Hypoxia',
  'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY' = 'ROS_Pathway',
  'HALLMARK_PI3K_AKT_MTOR_SIGNALING' = 'PI3K_AKT_mTOR'
)

# Get available genes in your Seurat object
available_genes <- rownames(so_subset)
cat("Total genes in dataset:", length(available_genes), "\n")

# Function to filter gene sets to only include available genes
filter_geneset <- function(geneset, available_genes, min_genes = 5) {
  overlap_genes <- intersect(geneset, available_genes)
  cat("Original:", length(geneset), "genes | Found:", length(overlap_genes), "genes")
  
  if(length(overlap_genes) >= min_genes) {
    cat(" ✓ KEEP\n")
    return(overlap_genes)
  } else {
    cat(" ✗ SKIP (too few genes)\n") 
    return(NULL)
  }
}

# Filter all HALLMARK gene sets
filtered_hallmark_genes <- list()
cat("\nFiltering HALLMARK pathways:\n")

for(pathway in kidney_hallmark_terms) {
  cat(pathway, ": ")
  filtered_genes <- filter_geneset(kidney_hallmark_genes[[pathway]], available_genes)
  
  if(!is.null(filtered_genes)) {
    filtered_hallmark_genes[[pathway]] <- filtered_genes
  }
}

# Now add module scores using filtered gene sets
cat("\nAdding module scores:\n")

for(pathway in names(filtered_hallmark_genes)) {
  score_name <- hallmark_names[pathway]
  
  so_subset <- AddModuleScore(
    object = so_subset,
    features = list(filtered_hallmark_genes[[pathway]]),
    name = score_name
  )
  
  cat("✓", score_name, ":", length(filtered_hallmark_genes[[pathway]]), "genes\n")
}
# Optional: Keep GO term analysis
go_term_names <- c(
  'GO0032868' = 'Response_to_Insulin',
  'GO0008286' = 'Insulin_Receptor_Signaling', 
  'GO0046627' = 'Neg_Reg_Insulin_Signaling',
  'GO0046628' = 'Pos_Reg_Insulin_Signaling',
  'GO0005159' = 'IGF_Receptor_Binding'
)

# Loop through and add each GO module score
for(i in 1:length(insulin_sensitivity_genes)) {
  go_id <- names(insulin_sensitivity_genes)[i]
  score_name <- go_term_names[str_replace(go_id, pattern = ':', replacement = '')]
  
  so_subset <- AddModuleScore(
    object = so_subset,
    features = list(insulin_sensitivity_genes[[i]]),
    name = score_name
  )
  
  cat("Added GO module score for", score_name, "with", length(insulin_sensitivity_genes[[i]]), "genes\n")
}




meta.data <- so_subset@meta.data
meta.data <- meta.data %>% 
  dplyr::select(record_id, mrn, group, celltype2, KPMP_celltype, epic_sglti2_1, 
                # Original custom scores
                TCA_score1, OxPhos_score1, 
                # GO terms (optional - if you kept them)
                Response_to_Insulin1, Insulin_Receptor_Signaling1, Neg_Reg_Insulin_Signaling1, 
                Pos_Reg_Insulin_Signaling1, IGF_Receptor_Binding1,
                # HALLMARK pathway scores - Metabolism
                OxPhos_Hallmark1, Fatty_Acid_Metabolism1, Glycolysis1, mTORC1_Signaling1, 
                Adipogenesis1, Peroxisome1,
                # HALLMARK pathway scores - Immune/Inflammatory  
                Inflammatory_Response1, IFN_Gamma_Response1, IFN_Alpha_Response1, 
                IL6_JAK_STAT31, TNFa_NF_kB1, Complement1, Allograft_Rejection1,
                # HALLMARK pathway scores - Cross-cutting
                Hypoxia1, ROS_Pathway1, PI3K_AKT_mTOR1)

write.table(meta.data, 'C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt', 
            row.names=F, quote=F, sep='\t')
remove(list=ls())



























##############module score graphing 
module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt')
module_scores$celltype2 <- NULL
module_scores$KPMP_celltype <- NULL





module_scores_summary <- module_scores %>%
  group_by(record_id, mrn, group) %>%
  summarize(
    # Calculate mean for all numeric module score columns
    across(c(TCA_score1, OxPhos_score1, 
             # GO terms
             Response_to_Insulin1, Insulin_Receptor_Signaling1, 
             Neg_Reg_Insulin_Signaling1, Pos_Reg_Insulin_Signaling1, 
             IGF_Receptor_Binding1,
             # HALLMARK pathways - Metabolism
             OxPhos_Hallmark1, Fatty_Acid_Metabolism1, Glycolysis1, 
             mTORC1_Signaling1, Adipogenesis1, Peroxisome1,
             # HALLMARK pathways - Immune/Inflammatory
             Inflammatory_Response1, IFN_Gamma_Response1, IFN_Alpha_Response1,
             IL6_JAK_STAT31, TNFa_NF_kB1, Complement1, Allograft_Rejection1,
             # HALLMARK pathways - Cross-cutting
             Hypoxia1, ROS_Pathway1, PI3K_AKT_mTOR1), 
           list(mean = ~mean(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = 'drop'
  )



combined_df <- module_scores_summary


########## Creating panel tables with new HALLMARK pathways

library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)
library(broom)

# Function to perform statistical tests and add significance annotations
add_significance <- function(data, variable, group_col = "group") {
  # Filter out NA values in both variable and group columns
  clean_data <- data %>% filter(!is.na(.data[[group_col]]) & !is.na(.data[[variable]]))
  
  # Check if we have enough data points
  if(nrow(clean_data) < 3) {
    return(list(significant = FALSE, p_value = NA))
  }
  
  # Get unique groups
  groups <- unique(clean_data[[group_col]])
  
  # If only one group, can't do statistical test
  if(length(groups) < 2) {
    return(list(significant = FALSE, p_value = NA))
  }
  
  # If exactly two groups, perform t-test
  if(length(groups) == 2) {
    group1_data <- clean_data[clean_data[[group_col]] == groups[1], variable]
    group2_data <- clean_data[clean_data[[group_col]] == groups[2], variable]
    
    # Perform unpaired t-test (assuming unequal variances)
    t_result <- t.test(group1_data, group2_data, var.equal = FALSE)
    
    return(list(significant = t_result$p.value < 0.05, p_value = t_result$p.value))
  }
  
  # If more than two groups, perform ANOVA (uncorrected)
  if(length(groups) > 2) {
    aov_result <- aov(as.formula(paste(variable, "~", group_col)), data = clean_data)
    p_value <- summary(aov_result)[[1]][["Pr(>F)"]][1]
    
    return(list(significant = p_value < 0.05, p_value = p_value))
  }
}

# Function to get pathway classification and color
get_pathway_color <- function(pathway_name) {
  custom_scores <- c("TCA_score", "OxPhos_score")
  go_terms <- c("Response_to_Insulin", "Insulin_Receptor_Signaling", 
                "Neg_Reg_Insulin_Signaling", "Pos_Reg_Insulin_Signaling", 
                "IGF_Receptor_Binding")
  metabolism_hallmark <- c("OxPhos_Hallmark", "Fatty_Acid_Metabolism", "Glycolysis", 
                           "mTORC1_Signaling", "Adipogenesis", "Peroxisome")
  immune_inflammatory <- c("Inflammatory_Response", "IFN_Gamma_Response", "IFN_Alpha_Response",
                           "IL6_JAK_STAT3", "TNFa_NF_kB", "Complement", "Allograft_Rejection")
  cross_cutting <- c("Hypoxia", "ROS_Pathway", "PI3K_AKT_mTOR")
  
  # Define colors for each category
  if(pathway_name %in% custom_scores) {
    return("#2E8B57")  # Sea Green for custom
  } else if(pathway_name %in% go_terms) {
    return("#4169E1")  # Royal Blue for GO terms
  } else if(pathway_name %in% metabolism_hallmark) {
    return("#FF6347")  # Tomato Red for metabolism
  } else if(pathway_name %in% immune_inflammatory) {
    return("#8A2BE2")  # Blue Violet for immune/inflammatory
  } else if(pathway_name %in% cross_cutting) {
    return("#FF8C00")  # Dark Orange for cross-cutting
  } else {
    return("#000000")  # Black for unknown
  }
}

# Function to create a single boxplot with significance
create_boxplot <- function(data, variable, title = NULL) {
  # Filter out NA values for plotting
  plot_data <- data %>% filter(!is.na(group) & !is.na(.data[[variable]]))
  
  # Check if we have enough data
  if(nrow(plot_data) < 3) {
    return(ggplot() + 
             theme_void() + 
             labs(title = paste("Insufficient data:", title)) +
             theme(plot.title = element_text(size = 10, hjust = 0.5, color = "gray")))
  }
  
  # Get significance results
  sig_results <- add_significance(data, variable)
  
  # Clean up title and get pathway classification
  clean_title <- if(is.null(title)) gsub("1_mean", "", variable) else title
  pathway_name <- gsub("_", " ", gsub("1_mean", "", variable))
  
  # Extract pathway name for color classification (remove spaces and formatting)
  pathway_for_color <- gsub(" ", "_", clean_title)
  pathway_for_color <- gsub("_Hallmark", "_Hallmark", pathway_for_color)  # Keep Hallmark suffix for classification
  if(pathway_for_color == "OxPhos") pathway_for_color <- "OxPhos_Hallmark"  # Handle the shortened version
  
  title_color <- get_pathway_color(pathway_for_color)
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = group, y = .data[[variable]], fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 10, hjust = 0.5, color = title_color),
      legend.position = "none"
    ) +
    labs(title = clean_title,
         y = "Score") +
    scale_fill_manual(values = c("Lean_Control" = "#619CFF", "Type_2_Diabetes" = "#F8766D"))
  
  # Add significance annotation if significant
  if(!is.na(sig_results$p_value) && sig_results$significant) {
    p <- p + annotate("text", x = Inf, y = Inf, 
                      label = paste0("p = ", round(sig_results$p_value, 4), "*"), 
                      hjust = 1.1, vjust = 1.1, size = 3, color = "red")
  }
  
  return(p)
}

# Function to create a legend plot
create_color_legend <- function() {
  legend_data <- data.frame(
    Category = c("Custom Scores", "GO Terms", "Metabolism", "Immune/Inflammatory", "Cross-cutting"),
    Color = c("#2E8B57", "#4169E1", "#FF6347", "#8A2BE2", "#FF8C00"),
    x = 1,
    y = c(5, 4, 3, 2, 1)
  )
  
  legend_plot <- ggplot(legend_data, aes(x = x, y = y)) +
    geom_text(aes(label = Category, color = Color), 
              size = 4, hjust = 0, fontface = "bold") +
    scale_color_identity() +
    xlim(0.5, 3) +
    ylim(0.5, 5.5) +
    theme_void() +
    labs(title = "Pathway Categories") +
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
  
  return(legend_plot)
}

# Function to create pathway-specific layouts
create_pathway_layout <- function(data, pathway_group = "all", include_legend = TRUE) {
  
  # Define pathway groups
  custom_scores <- c("TCA_score", "OxPhos_score")
  
  go_terms <- c("Response_to_Insulin", "Insulin_Receptor_Signaling", 
                "Neg_Reg_Insulin_Signaling", "Pos_Reg_Insulin_Signaling", 
                "IGF_Receptor_Binding")
  
  metabolism_hallmark <- c("OxPhos_Hallmark", "Fatty_Acid_Metabolism", "Glycolysis", 
                           "mTORC1_Signaling", "Adipogenesis", "Peroxisome")
  
  immune_inflammatory <- c("Inflammatory_Response", "IFN_Gamma_Response", "IFN_Alpha_Response",
                           "IL6_JAK_STAT3", "TNFa_NF_kB", "Complement", "Allograft_Rejection")
  
  cross_cutting <- c("Hypoxia", "ROS_Pathway", "PI3K_AKT_mTOR")
  
  # Select traits based on pathway group
  if(pathway_group == "all") {
    traits <- c(custom_scores, go_terms, metabolism_hallmark, immune_inflammatory, cross_cutting)
    plot_title <- "All Gene Set Enrichment Scores by Group - MEAN"
    ncols <- 5  # 5x5 grid for all pathways (23 pathways + 2 empty slots = 25)
    nrows <- 5
  } else if(pathway_group == "custom") {
    traits <- custom_scores
    plot_title <- "Custom Gene Set Scores - MEAN"
    ncols <- 2
    nrows <- 1
  } else if(pathway_group == "go_terms") {
    traits <- go_terms
    plot_title <- "GO Term Enrichment Scores - MEAN"
    ncols <- 3
    nrows <- 2
  } else if(pathway_group == "metabolism") {
    traits <- metabolism_hallmark
    plot_title <- "Metabolism HALLMARK Pathways - MEAN"
    ncols <- 3
    nrows <- 2
  } else if(pathway_group == "immune") {
    traits <- immune_inflammatory
    plot_title <- "Immune/Inflammatory HALLMARK Pathways - MEAN"
    ncols <- 3
    nrows <- 3
  } else if(pathway_group == "cross_cutting") {
    traits <- cross_cutting
    plot_title <- "Cross-cutting HALLMARK Pathways - MEAN"
    ncols <- 3
    nrows <- 1
  }
  
  # Create column names (assuming mean statistics)
  columns <- paste0(traits, "1_mean")
  
  # Create list of plots
  plots <- vector("list", length(traits))
  
  for(i in 1:length(traits)) {
    if(columns[i] %in% names(data)) {
      # Clean up title by removing underscores and making more readable
      clean_title <- gsub("_", " ", traits[i])
      clean_title <- gsub("Hallmark$", "", clean_title)  # Remove "Hallmark" suffix if present
      
      plots[[i]] <- create_boxplot(data, columns[i], title = clean_title)
    } else {
      # Create an empty placeholder plot for missing variables
      plots[[i]] <- ggplot() + 
        theme_void() + 
        labs(title = paste("Missing:", gsub("_", " ", traits[i]))) +
        theme(plot.title = element_text(size = 10, hjust = 0.5, color = "gray"))
    }
  }
  
  # Remove any NULL elements
  plots <- plots[!sapply(plots, is.null)]
  
  # Fill remaining slots for even grid if needed
  total_slots <- ncols * nrows
  if(length(plots) < total_slots) {
    for(i in (length(plots)+1):total_slots) {
      plots[[i]] <- ggplot() + theme_void()
    }
  }
  
  # Add legend if requested and if showing all pathways
  if(include_legend && pathway_group == "all") {
    legend_plot <- create_color_legend()
    
    # Create the main plots grid first
    main_plot <- wrap_plots(plots, ncol = ncols, nrow = nrows)
    
    # Combine main plot with legend using different approach
    final_plot <- main_plot | legend_plot
    final_plot <- final_plot + 
      plot_layout(widths = c(4, 1)) +
      plot_annotation(title = plot_title,
                      theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
  } else {
    # Combine into grid without legend
    combined_plot <- wrap_plots(plots, ncol = ncols, nrow = nrows)
    
    # Just add title without legend
    final_plot <- combined_plot + 
      plot_annotation(title = plot_title,
                      theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
  }
  
  return(final_plot)
}



all_pathways_plot <- create_pathway_layout(combined_df, 'all')


ggsave("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/allpathways_panel_comparisons.pdf", 
       all_pathways_plot, width = 25, height = 25, dpi = 300)


png("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/allpathways_panel_comparisons.png", 
       width = 25, height = 25, units = 'in', res = 300)
print(all_pathways_plot)
dev.off()



###############################################################################Graph again but with only PT cells 

#module score graphing 
module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt')
module_scores <- module_scores %>% 
  filter(celltype2 == 'PT')
module_scores$celltype2 <- NULL
module_scores$KPMP_celltype <- NULL




module_scores_summary <- module_scores %>%
  group_by(record_id, mrn, group) %>%
  summarize(
    # Calculate mean for all numeric module score columns
    across(c(TCA_score1, OxPhos_score1, 
             # GO terms
             Response_to_Insulin1, Insulin_Receptor_Signaling1, 
             Neg_Reg_Insulin_Signaling1, Pos_Reg_Insulin_Signaling1, 
             IGF_Receptor_Binding1,
             # HALLMARK pathways - Metabolism
             OxPhos_Hallmark1, Fatty_Acid_Metabolism1, Glycolysis1, 
             mTORC1_Signaling1, Adipogenesis1, Peroxisome1,
             # HALLMARK pathways - Immune/Inflammatory
             Inflammatory_Response1, IFN_Gamma_Response1, IFN_Alpha_Response1,
             IL6_JAK_STAT31, TNFa_NF_kB1, Complement1, Allograft_Rejection1,
             # HALLMARK pathways - Cross-cutting
             Hypoxia1, ROS_Pathway1, PI3K_AKT_mTOR1), 
           list(mean = ~mean(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = 'drop'
  )



combined_df <- module_scores_summary

#combined_df <- combined_df %>% filter(!is.na(avg_c_k2))

combined_df$mrn <- NULL


########## Creating panel tables with new HALLMARK pathways

library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)
library(broom)

# Function to perform statistical tests and add significance annotations
add_significance <- function(data, variable, group_col = "group") {
  # Filter out NA values in both variable and group columns
  clean_data <- data %>% filter(!is.na(.data[[group_col]]) & !is.na(.data[[variable]]))
  
  # Check if we have enough data points
  if(nrow(clean_data) < 3) {
    return(list(significant = FALSE, p_value = NA))
  }
  
  # Get unique groups
  groups <- unique(clean_data[[group_col]])
  
  # If only one group, can't do statistical test
  if(length(groups) < 2) {
    return(list(significant = FALSE, p_value = NA))
  }
  
  # If exactly two groups, perform t-test
  if(length(groups) == 2) {
    group1_data <- clean_data[clean_data[[group_col]] == groups[1], variable]
    group2_data <- clean_data[clean_data[[group_col]] == groups[2], variable]
    
    # Perform unpaired t-test (assuming unequal variances)
    t_result <- t.test(group1_data, group2_data, var.equal = FALSE)
    
    return(list(significant = t_result$p.value < 0.05, p_value = t_result$p.value))
  }
  
  # If more than two groups, perform ANOVA (uncorrected)
  if(length(groups) > 2) {
    aov_result <- aov(as.formula(paste(variable, "~", group_col)), data = clean_data)
    p_value <- summary(aov_result)[[1]][["Pr(>F)"]][1]
    
    return(list(significant = p_value < 0.05, p_value = p_value))
  }
}

# Function to get pathway classification and color
get_pathway_color <- function(pathway_name) {
  custom_scores <- c("TCA_score", "OxPhos_score")
  go_terms <- c("Response_to_Insulin", "Insulin_Receptor_Signaling", 
                "Neg_Reg_Insulin_Signaling", "Pos_Reg_Insulin_Signaling", 
                "IGF_Receptor_Binding")
  metabolism_hallmark <- c("OxPhos_Hallmark", "Fatty_Acid_Metabolism", "Glycolysis", 
                           "mTORC1_Signaling", "Adipogenesis", "Peroxisome")
  immune_inflammatory <- c("Inflammatory_Response", "IFN_Gamma_Response", "IFN_Alpha_Response",
                           "IL6_JAK_STAT3", "TNFa_NF_kB", "Complement", "Allograft_Rejection")
  cross_cutting <- c("Hypoxia", "ROS_Pathway", "PI3K_AKT_mTOR")
  
  # Define colors for each category
  if(pathway_name %in% custom_scores) {
    return("#2E8B57")  # Sea Green for custom
  } else if(pathway_name %in% go_terms) {
    return("#4169E1")  # Royal Blue for GO terms
  } else if(pathway_name %in% metabolism_hallmark) {
    return("#FF6347")  # Tomato Red for metabolism
  } else if(pathway_name %in% immune_inflammatory) {
    return("#8A2BE2")  # Blue Violet for immune/inflammatory
  } else if(pathway_name %in% cross_cutting) {
    return("#FF8C00")  # Dark Orange for cross-cutting
  } else {
    return("#000000")  # Black for unknown
  }
}

# Function to create a single boxplot with significance
create_boxplot <- function(data, variable, title = NULL) {
  # Filter out NA values for plotting
  plot_data <- data %>% filter(!is.na(group) & !is.na(.data[[variable]]))
  
  # Check if we have enough data
  if(nrow(plot_data) < 3) {
    return(ggplot() + 
             theme_void() + 
             labs(title = paste("Insufficient data:", title)) +
             theme(plot.title = element_text(size = 10, hjust = 0.5, color = "gray")))
  }
  
  # Get significance results
  sig_results <- add_significance(data, variable)
  
  # Clean up title and get pathway classification
  clean_title <- if(is.null(title)) gsub("1_mean", "", variable) else title
  pathway_name <- gsub("_", " ", gsub("1_mean", "", variable))
  
  # Extract pathway name for color classification (remove spaces and formatting)
  pathway_for_color <- gsub(" ", "_", clean_title)
  pathway_for_color <- gsub("_Hallmark", "_Hallmark", pathway_for_color)  # Keep Hallmark suffix for classification
  if(pathway_for_color == "OxPhos") pathway_for_color <- "OxPhos_Hallmark"  # Handle the shortened version
  
  title_color <- get_pathway_color(pathway_for_color)
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = group, y = .data[[variable]], fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 10, hjust = 0.5, color = title_color),
      legend.position = "none"
    ) +
    labs(title = clean_title,
         y = "Score") +
    scale_fill_manual(values = c("Lean_Control" = "#619CFF", "Type_2_Diabetes" = "#F8766D"))
  
  # Add significance annotation if significant
  if(!is.na(sig_results$p_value) && sig_results$significant) {
    p <- p + annotate("text", x = Inf, y = Inf, 
                      label = paste0("p = ", round(sig_results$p_value, 4), "*"), 
                      hjust = 1.1, vjust = 1.1, size = 3, color = "red")
  }
  
  return(p)
}

# Function to create a legend plot
create_color_legend <- function() {
  legend_data <- data.frame(
    Category = c("Custom Scores", "GO Terms", "Metabolism", "Immune/Inflammatory", "Cross-cutting"),
    Color = c("#2E8B57", "#4169E1", "#FF6347", "#8A2BE2", "#FF8C00"),
    x = 1,
    y = c(5, 4, 3, 2, 1)
  )
  
  legend_plot <- ggplot(legend_data, aes(x = x, y = y)) +
    geom_text(aes(label = Category, color = Color), 
              size = 4, hjust = 0, fontface = "bold") +
    scale_color_identity() +
    xlim(0.5, 3) +
    ylim(0.5, 5.5) +
    theme_void() +
    labs(title = "Pathway Categories") +
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
  
  return(legend_plot)
}

# Function to create pathway-specific layouts
create_pathway_layout <- function(data, pathway_group = "all", include_legend = TRUE) {
  
  # Define pathway groups
  custom_scores <- c("TCA_score", "OxPhos_score")
  
  go_terms <- c("Response_to_Insulin", "Insulin_Receptor_Signaling", 
                "Neg_Reg_Insulin_Signaling", "Pos_Reg_Insulin_Signaling", 
                "IGF_Receptor_Binding")
  
  metabolism_hallmark <- c("OxPhos_Hallmark", "Fatty_Acid_Metabolism", "Glycolysis", 
                           "mTORC1_Signaling", "Adipogenesis", "Peroxisome")
  
  immune_inflammatory <- c("Inflammatory_Response", "IFN_Gamma_Response", "IFN_Alpha_Response",
                           "IL6_JAK_STAT3", "TNFa_NF_kB", "Complement", "Allograft_Rejection")
  
  cross_cutting <- c("Hypoxia", "ROS_Pathway", "PI3K_AKT_mTOR")
  
  # Select traits based on pathway group
  if(pathway_group == "all") {
    traits <- c(custom_scores, go_terms, metabolism_hallmark, immune_inflammatory, cross_cutting)
    plot_title <- "All Gene Set Enrichment Scores by Group - MEAN"
    ncols <- 5  # 5x5 grid for all pathways (23 pathways + 2 empty slots = 25)
    nrows <- 5
  } else if(pathway_group == "custom") {
    traits <- custom_scores
    plot_title <- "Custom Gene Set Scores - MEAN"
    ncols <- 2
    nrows <- 1
  } else if(pathway_group == "go_terms") {
    traits <- go_terms
    plot_title <- "GO Term Enrichment Scores - MEAN"
    ncols <- 3
    nrows <- 2
  } else if(pathway_group == "metabolism") {
    traits <- metabolism_hallmark
    plot_title <- "Metabolism HALLMARK Pathways - MEAN"
    ncols <- 3
    nrows <- 2
  } else if(pathway_group == "immune") {
    traits <- immune_inflammatory
    plot_title <- "Immune/Inflammatory HALLMARK Pathways - MEAN"
    ncols <- 3
    nrows <- 3
  } else if(pathway_group == "cross_cutting") {
    traits <- cross_cutting
    plot_title <- "Cross-cutting HALLMARK Pathways - MEAN"
    ncols <- 3
    nrows <- 1
  }
  
  # Create column names (assuming mean statistics)
  columns <- paste0(traits, "1_mean")
  
  # Create list of plots
  plots <- vector("list", length(traits))
  
  for(i in 1:length(traits)) {
    if(columns[i] %in% names(data)) {
      # Clean up title by removing underscores and making more readable
      clean_title <- gsub("_", " ", traits[i])
      clean_title <- gsub("Hallmark$", "", clean_title)  # Remove "Hallmark" suffix if present
      
      plots[[i]] <- create_boxplot(data, columns[i], title = clean_title)
    } else {
      # Create an empty placeholder plot for missing variables
      plots[[i]] <- ggplot() + 
        theme_void() + 
        labs(title = paste("Missing:", gsub("_", " ", traits[i]))) +
        theme(plot.title = element_text(size = 10, hjust = 0.5, color = "gray"))
    }
  }
  
  # Remove any NULL elements
  plots <- plots[!sapply(plots, is.null)]
  
  # Fill remaining slots for even grid if needed
  total_slots <- ncols * nrows
  if(length(plots) < total_slots) {
    for(i in (length(plots)+1):total_slots) {
      plots[[i]] <- ggplot() + theme_void()
    }
  }
  
  # Add legend if requested and if showing all pathways
  if(include_legend && pathway_group == "all") {
    legend_plot <- create_color_legend()
    
    # Create the main plots grid first
    main_plot <- wrap_plots(plots, ncol = ncols, nrow = nrows)
    
    # Combine main plot with legend using different approach
    final_plot <- main_plot | legend_plot
    final_plot <- final_plot + 
      plot_layout(widths = c(4, 1)) +
      plot_annotation(title = plot_title,
                      theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
  } else {
    # Combine into grid without legend
    combined_plot <- wrap_plots(plots, ncol = ncols, nrow = nrows)
    
    # Just add title without legend
    final_plot <- combined_plot + 
      plot_annotation(title = plot_title,
                      theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
  }
  
  return(final_plot)
}



all_pathways_plot <- create_pathway_layout(combined_df, 'all')


ggsave("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/allpathways_panel_comparisons_PTCells.pdf", 
       all_pathways_plot, width = 25, height = 25, dpi = 300)





























############################################## Gene module scores with PET variables 




#module score graphing 
module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt')
module_scores <- module_scores %>% 
  filter(celltype2 == 'PT')
module_scores$celltype2 <- NULL
module_scores$KPMP_celltype <- NULL




module_scores_summary <- module_scores %>%
  group_by(record_id, mrn, group) %>%
  summarize(
    # Calculate mean for all numeric module score columns
    across(c(TCA_score1, OxPhos_score1, 
             # GO terms
             Response_to_Insulin1, Insulin_Receptor_Signaling1, 
             Neg_Reg_Insulin_Signaling1, Pos_Reg_Insulin_Signaling1, 
             IGF_Receptor_Binding1,
             # HALLMARK pathways - Metabolism
             OxPhos_Hallmark1, Fatty_Acid_Metabolism1, Glycolysis1, 
             mTORC1_Signaling1, Adipogenesis1, Peroxisome1,
             # HALLMARK pathways - Immune/Inflammatory
             Inflammatory_Response1, IFN_Gamma_Response1, IFN_Alpha_Response1,
             IL6_JAK_STAT31, TNFa_NF_kB1, Complement1, Allograft_Rejection1,
             # HALLMARK pathways - Cross-cutting
             Hypoxia1, ROS_Pathway1, PI3K_AKT_mTOR1), 
           list(mean = ~mean(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = 'drop'
  )



combined_df <- module_scores_summary





#label(combined_df$hba1c) <- "HbA1c (%)"
#label(combined_df$age) <- "Age (Years)"
#label(combined_df$sex) <- "Sex"
#label(combined_df$race_ethnicity) <- "Race/Ethnicity"
#label(combined_df$bmi) <- "BMI (kg/m2)"
#label(combined_df$epic_sglti2_1) <- 'SGLT2 Inhibitor Use'

#table1::table1(~age + sex + race_ethnicity + bmi +hba1c +  study + epic_sglti2_1| group, 
#               data = combined_df)


# Fix data types before creating the table
library(gtsummary)
library(gt)
library(dplyr)

# Convert variables to proper data types
combined_df <- combined_df %>%
  mutate(
    # Ensure continuous variables are numeric
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    hba1c = as.numeric(hba1c),
    
    # Ensure categorical variables are factors or characters
    sex = as.factor(sex),
    race_ethnicity = as.factor(race_ethnicity),
    study = as.factor(study),
    group = as.factor(group),
    epic_sglti2_1 = as.factor(epic_sglti2_1)
  )



# Now create the table with proper data types
desc_table1_fixed <- combined_df %>%
  select(age, sex, race_ethnicity, bmi, hba1c, study, group, epic_sglti2_1) %>%
  tbl_summary(
    by = group,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      hba1c ~ "continuous",
      sex ~ "categorical",
      race_ethnicity ~ "categorical",
      study ~ "categorical",
      epic_sglti2_1 ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      sex ~ "Sex", 
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      study ~ "Study",
      epic_sglti2_1 ~ "SGLT2 Inhibitor Use"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "t.test"
    # Skip categorical p-values if they cause issues
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Group**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

# Save version with epic
desc_table1_fixed %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/demographics_aim5_with_epic_final.png", 
         vwidth = 1200, vheight = 800)



#plotting
analysis_df <- combined_df %>% 
  dplyr::select(2:49, 51:63)

# Function to create correlation matrix and plot (CORRECTED VERSION)
create_correlation_heatmap <- function(data, pet_variables, module_variables, stat_type, subset_name = "all", save_plot = TRUE) {
  
  # Select the variables we need
  selected_vars <- c(module_variables, pet_variables)
  correlation_data <- data[, selected_vars]
  
  # Calculate correlation between module scores (rows) and PET variables (columns)
  # We need to transpose because we want modules as rows and PET as columns
  corr_matrix <- cor(correlation_data[, module_variables], 
                     correlation_data[, pet_variables], 
                     use = 'pairwise.complete.obs')
  
  # Calculate p-values for the rectangular correlation matrix
  n_modules <- length(module_variables)
  n_pet <- length(pet_variables)
  p_mat <- matrix(NA, nrow = n_modules, ncol = n_pet)
  
  for (i in 1:n_modules) {
    for (j in 1:n_pet) {
      valid_pairs <- complete.cases(correlation_data[, module_variables[i]], 
                                    correlation_data[, pet_variables[j]])
      if (sum(valid_pairs) >= 3) {
        tryCatch({
          tmp <- cor.test(correlation_data[, module_variables[i]], 
                          correlation_data[, pet_variables[j]])
          p_mat[i, j] <- tmp$p.value
        }, error = function(e) {
          p_mat[i, j] <<- NA
        })
      }
    }
  }
  
  # Set row and column names
  rownames(corr_matrix) <- module_variables
  colnames(corr_matrix) <- pet_variables
  rownames(p_mat) <- module_variables
  colnames(p_mat) <- pet_variables
  
  # Clean up row names (remove _mean, _median, _sd suffixes)
  clean_row_names <- gsub("_mean$|_median$|_sd$", "", rownames(corr_matrix))
  clean_row_names <- gsub("1$", "", clean_row_names)  # Remove trailing "1"
  clean_row_names <- gsub("_", " ", clean_row_names)  # Replace underscores with spaces
  
  # Clean up column names for PET variables with specific labels
  pet_labels <- c(
    "avg_c_f" = "Cortical F",
    "avg_c_k2" = "Cortical K2", 
    "avg_c_k2_f" = "Cortical K2/F",
    "avg_c_f_vw" = "Cortical F (voxel)",
    "avg_c_k2_vw" = "Cortical K2 (voxel)",
    "avg_c_k2_f_vw" = "Cortical K2/F (voxel)"
  )
  
  clean_col_names <- pet_labels[colnames(corr_matrix)]
  
  rownames(corr_matrix) <- clean_row_names
  colnames(corr_matrix) <- clean_col_names
  rownames(p_mat) <- clean_row_names
  colnames(p_mat) <- clean_col_names
  
  # Create the plot
  if (save_plot) {
    png(paste0("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/pet_module_correlation_", stat_type, "_", subset_name, ".png"), 
        width = 1200, height = 800, res = 150)
  }
  
  # Plot with significance indicators - using safer approach
  # Check if we have any significant correlations
  has_significant <- any(p_mat < 0.05, na.rm = TRUE)
  
  # Create title with subset information
  title_text <- paste("PET Scan vs Module Scores Correlation -", toupper(stat_type))
  if (subset_name == "t2d") {
    title_text <- paste(title_text, "(Type 2 Diabetes Only)")
  } else if (subset_name == "all") {
    title_text <- paste(title_text, "(All Participants)")
  }
  
  if (has_significant) {
    # Plot with significance indicators
    corrplot(corr_matrix, 
             method = "color",
             type = "full",
             tl.col = "black",
             tl.cex = 0.8,
             tl.srt = 45,
             cl.cex = 0.8,
             title = title_text,
             mar = c(0, 0, 2, 0),
             p.mat = p_mat,
             sig.level = 0.05,
             insig = "pch",  # Use pch instead of label_sig
             pch = 4,        # Use 'x' for non-significant
             pch.cex = 0.8,
             pch.col = "white")
  } else {
    # Plot without significance indicators if none are significant
    corrplot(corr_matrix, 
             method = "color",
             type = "full",
             tl.col = "black",
             tl.cex = 0.8,
             tl.srt = 45,
             cl.cex = 0.8,
             title = title_text,
             mar = c(0, 0, 2, 0))
  }
  
  if (save_plot) {
    dev.off()
    cat("Saved:", paste0("pet_module_correlation_", stat_type, "_", subset_name, ".png\n"))
  }
  
  # Return the matrices for further analysis if needed
  return(list(correlation = corr_matrix, p_values = p_mat))
}

# Check what values are in the group variable
print("Group variable values:")
print(table(analysis_df$group, useNA = "ifany"))

# Create datasets for analysis
all_data <- analysis_df
t2d_data <- analysis_df[analysis_df$group == "Type 2 Diabetes", ]

print(paste("All participants: N =", nrow(all_data)))
print(paste("T2D participants: N =", nrow(t2d_data)))

# Define PET scan variables (columns)
pet_vars <- c("avg_c_f", "avg_c_k2", "avg_c_k2_f", "avg_c_f_vw", "avg_c_k2_vw", "avg_c_k2_f_vw")

# Get module score variables for each type
module_vars <- list(
  mean = grep("_mean$", names(analysis_df), value = TRUE),
  median = grep("_median$", names(analysis_df), value = TRUE),
  sd = grep("_sd$", names(analysis_df), value = TRUE)
)

# Remove "group" from module variables if it exists
module_vars <- lapply(module_vars, function(x) x[!grepl("group", x)])

print("Creating correlation heatmaps...")

# === ALL PARTICIPANTS ===
print("\n=== CREATING HEATMAPS FOR ALL PARTICIPANTS ===")

# 1. Mean correlations - All participants
print("Creating MEAN correlation heatmap (All participants)...")
mean_results_all <- create_correlation_heatmap(all_data, pet_vars, module_vars$mean, "mean", "all")

# 2. Median correlations - All participants
print("Creating MEDIAN correlation heatmap (All participants)...")
median_results_all <- create_correlation_heatmap(all_data, pet_vars, module_vars$median, "median", "all")

# 3. Standard deviation correlations - All participants
print("Creating SD correlation heatmap (All participants)...")
sd_results_all <- create_correlation_heatmap(all_data, pet_vars, module_vars$sd, "sd", "all")

# === TYPE 2 DIABETES PARTICIPANTS ONLY ===
print("\n=== CREATING HEATMAPS FOR T2D PARTICIPANTS ONLY ===")

# 1. Mean correlations - T2D only
print("Creating MEAN correlation heatmap (T2D only)...")
mean_results_t2d <- create_correlation_heatmap(t2d_data, pet_vars, module_vars$mean, "mean", "t2d")

# 2. Median correlations - T2D only
print("Creating MEDIAN correlation heatmap (T2D only)...")
median_results_t2d <- create_correlation_heatmap(t2d_data, pet_vars, module_vars$median, "median", "t2d")

# 3. Standard deviation correlations - T2D only
print("Creating SD correlation heatmap (T2D only)...")
sd_results_t2d <- create_correlation_heatmap(t2d_data, pet_vars, module_vars$sd, "sd", "t2d")

print("\n=== ALL HEATMAPS CREATED SUCCESSFULLY! ===")

# Summary of files created
cat("\nFiles created:\n")
cat("ALL PARTICIPANTS:\n")
cat("  - pet_module_correlation_mean_all.png\n")
cat("  - pet_module_correlation_median_all.png\n")
cat("  - pet_module_correlation_sd_all.png\n")
cat("\nTYPE 2 DIABETES ONLY:\n")
cat("  - pet_module_correlation_mean_t2d.png\n")
cat("  - pet_module_correlation_median_t2d.png\n")
cat("  - pet_module_correlation_sd_t2d.png\n")

# Optional: Print comparison summary of significant correlations
print("\n=== COMPARISON OF SIGNIFICANT CORRELATIONS ===")

# Function to count significant correlations
count_significant <- function(results) {
  sum(results$p_values < 0.05, na.rm = TRUE)
}

# Compare results
cat("\nNumber of significant correlations (p < 0.05):\n")
cat("MEAN:\n")
cat("  All participants:", count_significant(mean_results_all), "\n")
cat("  T2D only:", count_significant(mean_results_t2d), "\n")
cat("MEDIAN:\n")
cat("  All participants:", count_significant(median_results_all), "\n")
cat("  T2D only:", count_significant(median_results_t2d), "\n")
cat("SD:\n")
cat("  All participants:", count_significant(sd_results_all), "\n")
cat("  T2D only:", count_significant(sd_results_t2d), "\n")



######################################################################################################################








### Only PT Cells 

module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/GO_pathways_modulescores.txt')
module_scores <- module_scores %>% filter(celltype2 == 'PT')
module_scores$celltype2 <- NULL
module_scores$KPMP_celltype <- NULL


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

dat_results$avg_c_k2 <- NULL
dat_results$avg_m_k2 <- NULL
dat_results$avg_c_f <- NULL


dat_results <- dat_results %>% bind_cols(tmp_results, tmp_results_vw)


dat_results <- dat_results %>% filter(!is.na(avg_c_k2))

dat_results <- dat_results %>% filter(group %in% c('Lean Control', 'Type 2 Diabetes'))


dat_results <- dat_results %>% 
  dplyr::select(mrn, group, starts_with('avg_c'), age, sex, bmi, hba1c, study, epic_sglti2_1, race_ethnicity)




module_scores_summary <- module_scores %>%
  group_by(mrn) %>%
  summarize(
    # Calculate mean, median, and SD for all numeric columns
    across(c(TCA_score1, OxPhos_score1, Response_to_Insulin1, 
             Insulin_Receptor_Signaling1, Neg_Reg_Insulin_Signaling1, 
             Pos_Reg_Insulin_Signaling1, IGF_Receptor_Binding1,
             Chloride_Homeostasis1, Potassium_Homeostasis1,
             Sodium_Homeostasis1, Anion_Homeostasis1,
             Pos_Reg_Phosphorus_Metabolism1, Pos_Reg_Phosphate_Metabolism1,
             Protein_Catabolic_Regulation1, Renal_Sodium_Transport1,
             Renal_Sodium_Absorption1), 
           list(mean = ~mean(.x, na.rm = TRUE),
                median = ~median(.x, na.rm = TRUE), 
                sd = ~sd(.x, na.rm = TRUE)), 
           .names = "{.col}_{.fn}"),
    .groups = 'drop'
  )


combined_df <- module_scores_summary %>% 
  left_join(dat_results, by=c('mrn'))

combined_df <- combined_df %>% filter(!is.na(avg_c_k2))



#label(combined_df$hba1c) <- "HbA1c (%)"
#label(combined_df$age) <- "Age (Years)"
#label(combined_df$sex) <- "Sex"
#label(combined_df$race_ethnicity) <- "Race/Ethnicity"
#label(combined_df$bmi) <- "BMI (kg/m2)"
#label(combined_df$epic_sglti2_1) <- 'SGLT2 Inhibitor Use'

#table1::table1(~age + sex + race_ethnicity + bmi +hba1c +  study + epic_sglti2_1| group, 
#               data = combined_df)


# Fix data types before creating the table
library(gtsummary)
library(gt)
library(dplyr)

# Convert variables to proper data types
combined_df <- combined_df %>%
  mutate(
    # Ensure continuous variables are numeric
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    hba1c = as.numeric(hba1c),
    
    # Ensure categorical variables are factors or characters
    sex = as.factor(sex),
    race_ethnicity = as.factor(race_ethnicity),
    study = as.factor(study),
    group = as.factor(group),
    epic_sglti2_1 = as.factor(epic_sglti2_1)
  )



#plotting
analysis_df <- combined_df %>% 
  dplyr::select(2:62)

# Function to create correlation matrix and plot (CORRECTED VERSION)
create_correlation_heatmap <- function(data, pet_variables, module_variables, stat_type, subset_name = "all", save_plot = TRUE) {
  
  # Select the variables we need
  selected_vars <- c(module_variables, pet_variables)
  correlation_data <- data[, selected_vars]
  
  # Calculate correlation between module scores (rows) and PET variables (columns)
  # We need to transpose because we want modules as rows and PET as columns
  corr_matrix <- cor(correlation_data[, module_variables], 
                     correlation_data[, pet_variables], 
                     use = 'pairwise.complete.obs')
  
  # Calculate p-values for the rectangular correlation matrix
  n_modules <- length(module_variables)
  n_pet <- length(pet_variables)
  p_mat <- matrix(NA, nrow = n_modules, ncol = n_pet)
  
  for (i in 1:n_modules) {
    for (j in 1:n_pet) {
      valid_pairs <- complete.cases(correlation_data[, module_variables[i]], 
                                    correlation_data[, pet_variables[j]])
      if (sum(valid_pairs) >= 3) {
        tryCatch({
          tmp <- cor.test(correlation_data[, module_variables[i]], 
                          correlation_data[, pet_variables[j]])
          p_mat[i, j] <- tmp$p.value
        }, error = function(e) {
          p_mat[i, j] <<- NA
        })
      }
    }
  }
  
  # Set row and column names
  rownames(corr_matrix) <- module_variables
  colnames(corr_matrix) <- pet_variables
  rownames(p_mat) <- module_variables
  colnames(p_mat) <- pet_variables
  
  # Clean up row names (remove _mean, _median, _sd suffixes)
  clean_row_names <- gsub("_mean$|_median$|_sd$", "", rownames(corr_matrix))
  clean_row_names <- gsub("1$", "", clean_row_names)  # Remove trailing "1"
  clean_row_names <- gsub("_", " ", clean_row_names)  # Replace underscores with spaces
  
  # Clean up column names for PET variables with specific labels
  pet_labels <- c(
    "avg_c_f" = "Cortical F",
    "avg_c_k2" = "Cortical K2", 
    "avg_c_k2_f" = "Cortical K2/F",
    "avg_c_f_vw" = "Cortical F (voxel)",
    "avg_c_k2_vw" = "Cortical K2 (voxel)",
    "avg_c_k2_f_vw" = "Cortical K2/F (voxel)"
  )
  
  clean_col_names <- pet_labels[colnames(corr_matrix)]
  
  rownames(corr_matrix) <- clean_row_names
  colnames(corr_matrix) <- clean_col_names
  rownames(p_mat) <- clean_row_names
  colnames(p_mat) <- clean_col_names
  
  # Create the plot
  if (save_plot) {
    png(paste0("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/PT_cell_pet_module_correlation_", stat_type, "_", subset_name, ".png"), 
        width = 1200, height = 800, res = 150)
  }
  
  # Plot with significance indicators - using safer approach
  # Check if we have any significant correlations
  has_significant <- any(p_mat < 0.05, na.rm = TRUE)
  
  # Create title with subset information
  title_text <- paste("PET Scan vs Module Scores Correlation in PT Cells -", toupper(stat_type))
  if (subset_name == "t2d") {
    title_text <- paste(title_text, "(Type 2 Diabetes Only)")
  } else if (subset_name == "all") {
    title_text <- paste(title_text, "(All Participants)")
  }
  
  if (has_significant) {
    # Plot with significance indicators
    corrplot(corr_matrix, 
             method = "color",
             type = "full",
             tl.col = "black",
             tl.cex = 0.8,
             tl.srt = 45,
             cl.cex = 0.8,
             title = title_text,
             mar = c(0, 0, 2, 0),
             p.mat = p_mat,
             sig.level = 0.05,
             insig = "pch",  # Use pch instead of label_sig
             pch = 4,        # Use 'x' for non-significant
             pch.cex = 0.8,
             pch.col = "white")
  } else {
    # Plot without significance indicators if none are significant
    corrplot(corr_matrix, 
             method = "color",
             type = "full",
             tl.col = "black",
             tl.cex = 0.8,
             tl.srt = 45,
             cl.cex = 0.8,
             title = title_text,
             mar = c(0, 0, 2, 0))
  }
  
  if (save_plot) {
    dev.off()
    cat("Saved:", paste0("pet_module_correlation_", stat_type, "_", subset_name, ".png\n"))
  }
  
  # Return the matrices for further analysis if needed
  return(list(correlation = corr_matrix, p_values = p_mat))
}

# Check what values are in the group variable
print("Group variable values:")
print(table(analysis_df$group, useNA = "ifany"))

# Create datasets for analysis
all_data <- analysis_df
t2d_data <- analysis_df[analysis_df$group == "Type 2 Diabetes", ]

print(paste("All participants: N =", nrow(all_data)))
print(paste("T2D participants: N =", nrow(t2d_data)))

# Define PET scan variables (columns)
pet_vars <- c("avg_c_f", "avg_c_k2", "avg_c_k2_f", "avg_c_f_vw", "avg_c_k2_vw", "avg_c_k2_f_vw")

# Get module score variables for each type
module_vars <- list(
  mean = grep("_mean$", names(analysis_df), value = TRUE),
  median = grep("_median$", names(analysis_df), value = TRUE),
  sd = grep("_sd$", names(analysis_df), value = TRUE)
)

# Remove "group" from module variables if it exists
module_vars <- lapply(module_vars, function(x) x[!grepl("group", x)])

print("Creating correlation heatmaps...")

# === ALL PARTICIPANTS ===
print("\n=== CREATING HEATMAPS FOR ALL PARTICIPANTS ===")

# 1. Mean correlations - All participants
print("Creating MEAN correlation heatmap (All participants)...")
mean_results_all <- create_correlation_heatmap(all_data, pet_vars, module_vars$mean, "mean", "all")

# 2. Median correlations - All participants
print("Creating MEDIAN correlation heatmap (All participants)...")
median_results_all <- create_correlation_heatmap(all_data, pet_vars, module_vars$median, "median", "all")

# 3. Standard deviation correlations - All participants
print("Creating SD correlation heatmap (All participants)...")
sd_results_all <- create_correlation_heatmap(all_data, pet_vars, module_vars$sd, "sd", "all")

# === TYPE 2 DIABETES PARTICIPANTS ONLY ===
print("\n=== CREATING HEATMAPS FOR T2D PARTICIPANTS ONLY ===")

# 1. Mean correlations - T2D only
print("Creating MEAN correlation heatmap (T2D only)...")
mean_results_t2d <- create_correlation_heatmap(t2d_data, pet_vars, module_vars$mean, "mean", "t2d")

# 2. Median correlations - T2D only
print("Creating MEDIAN correlation heatmap (T2D only)...")
median_results_t2d <- create_correlation_heatmap(t2d_data, pet_vars, module_vars$median, "median", "t2d")

# 3. Standard deviation correlations - T2D only
print("Creating SD correlation heatmap (T2D only)...")
sd_results_t2d <- create_correlation_heatmap(t2d_data, pet_vars, module_vars$sd, "sd", "t2d")

print("\n=== ALL HEATMAPS CREATED SUCCESSFULLY! ===")

# Summary of files created
cat("\nFiles created:\n")
cat("ALL PARTICIPANTS:\n")
cat("  - pet_module_correlation_mean_all.png\n")
cat("  - pet_module_correlation_median_all.png\n")
cat("  - pet_module_correlation_sd_all.png\n")
cat("\nTYPE 2 DIABETES ONLY:\n")
cat("  - pet_module_correlation_mean_t2d.png\n")
cat("  - pet_module_correlation_median_t2d.png\n")
cat("  - pet_module_correlation_sd_t2d.png\n")

# Optional: Print comparison summary of significant correlations
print("\n=== COMPARISON OF SIGNIFICANT CORRELATIONS ===")

# Function to count significant correlations
count_significant <- function(results) {
  sum(results$p_values < 0.05, na.rm = TRUE)
}

# Compare results
cat("\nNumber of significant correlations (p < 0.05):\n")
cat("MEAN:\n")
cat("  All participants:", count_significant(mean_results_all), "\n")
cat("  T2D only:", count_significant(mean_results_t2d), "\n")
cat("MEDIAN:\n")
cat("  All participants:", count_significant(median_results_all), "\n")
cat("  T2D only:", count_significant(median_results_t2d), "\n")
cat("SD:\n")
cat("  All participants:", count_significant(sd_results_all), "\n")
cat("  T2D only:", count_significant(sd_results_t2d), "\n")





















################ Plotting Correlations with PET variables 

remove(list=ls())

module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt')

celltype <- 'PT'



if(celltype == 'PT'){
module_scores <- module_scores %>% 
  filter(celltype2 == 'PT')
cell_text <- 'PT'
}else{
  cell_text <- 'All'
}
module_scores$celltype2 <- NULL
module_scores$KPMP_celltype <- NULL





module_scores_summary <- module_scores %>%
  group_by(record_id, mrn, group) %>%
  summarize(
    # Calculate mean for all numeric module score columns
    across(c(TCA_score1, OxPhos_score1, 
             # GO terms
             Response_to_Insulin1, Insulin_Receptor_Signaling1, 
             Neg_Reg_Insulin_Signaling1, Pos_Reg_Insulin_Signaling1, 
             IGF_Receptor_Binding1,
             # HALLMARK pathways - Metabolism
             OxPhos_Hallmark1, Fatty_Acid_Metabolism1, Glycolysis1, 
             mTORC1_Signaling1, Adipogenesis1, Peroxisome1,
             # HALLMARK pathways - Immune/Inflammatory
             Inflammatory_Response1, IFN_Gamma_Response1, IFN_Alpha_Response1,
             IL6_JAK_STAT31, TNFa_NF_kB1, Complement1, Allograft_Rejection1,
             # HALLMARK pathways - Cross-cutting
             Hypoxia1, ROS_Pathway1, PI3K_AKT_mTOR1), 
           list(mean = ~mean(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = 'drop'
  )




harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

#date_of_screen
#screen_date

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))

#gbm <- gbm %>% left_join(dat %>% dplyr::select(record_id, mrn), by='record_id')


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

dat_results$avg_c_k2 <- NULL
dat_results$avg_c_f <- NULL

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


# Function to get pathway colors (same as before)
get_pathway_color <- function(pathway_name) {
  custom_scores <- c("TCA", "OxPhos")
  go_terms <- c("Response to Insulin", "Insulin Receptor Signaling", 
                "Negative Regulation of Insulin Signaling", "Positive Regulation of Insulin Signaling", 
                "IGF Receptor Binding")
  metabolism_hallmark <- c("OxPhos (Hallmark)", "Fatty Acid Metabolism", "Glycolysis", 
                           "mTORC1 Signaling", "Adipogenesis", "Peroxisome")
  immune_inflammatory <- c("Inflammatory Response", "IFN Gamma Response", "IFN Alpha Response",
                           "IL6 JAK/STAT3", "TNFa NF kB1", "Complement", "Allograft Rejection")
  cross_cutting <- c("Hypoxia", "ROS Pathway", "PI3K AKT mTOR1")
  
  if(pathway_name %in% custom_scores) {
    return("#2E8B57")  # Sea Green
  } else if(pathway_name %in% go_terms) {
    return("#4169E1")  # Royal Blue
  } else if(pathway_name %in% metabolism_hallmark) {
    return("#FF6347")  # Tomato Red
  } else if(pathway_name %in% immune_inflammatory) {
    return("#8A2BE2")  # Blue Violet
  } else if(pathway_name %in% cross_cutting) {
    return("#FF8C00")  # Dark Orange
  } else {
    return("#000000")  # Black for unknown
  }
}

combined_df <- module_scores_summary %>% left_join(dat2 %>% 
                                                     dplyr::select(mrn, avg_c_k2, avg_c_f, avg_c_k2_f, 
                                                                   avg_c_k2_vw, avg_c_f_vw, avg_c_k2_f_vw), by= 'mrn')




# Filter to T2D only
t2d_only <- combined_df %>% 
  filter(group == 'Type_2_Diabetes') %>%
  filter(!is.na(avg_c_k2))

lc_only <- combined_df %>% 
  filter(group == 'Lean_Control') %>%
  filter(!is.na(avg_c_k2))

combined_df <- combined_df %>% dplyr::select(-record_id, -mrn, -group)


library(corrplot)

# Calculate correlations
combined_df_corr <- cor(combined_df, use = 'pairwise.complete.obs', method = 'spearman')

# Calculate p-values using cor.mtest
p_values <- cor.mtest(combined_df, method = 'spearman')

# Create subset for plotting
corr_subset <- as.matrix(combined_df_corr[c(1:23), c(24:29), drop = F])
p_subset <- as.matrix(p_values$p[c(1:23), c(24:29), drop = F])

# Add meaningful names
row_names <- c('TCA', 'OxPhos', 'Response to Insulin', 
               'Insulin Receptor Signaling', 'Negative Regulation of Insulin Signaling', 
               'Positive Regulation of Insulin Signaling', 'IGF Receptor Binding', 
               'OxPhos (Hallmark)', 'Fatty Acid Metabolism', 'Glycolysis', 'mTORC1 Signaling', 
               'Adipogenesis', 'Peroxisome', 'Inflammatory Response', 'IFN Gamma Response', 
               'IFN Alpha Response', 'IL6 JAK/STAT3', 'TNFa NF kB1', 'Complement', 
               'Allograft Rejection', 'Hypoxia', 'ROS Pathway', 'PI3K AKT mTOR1')

col_names <- c('Cortical K2', 'Cortical F', 'Cortical K2/F', 
               'Cortical K2 (voxel)', 'Cortical F (voxel)', 'Cortical K2/F (voxel)')

rownames(corr_subset) <- row_names
colnames(corr_subset) <- col_names

# Apply same names to p-value matrix
rownames(p_subset) <- rownames(corr_subset)
colnames(p_subset) <- colnames(corr_subset)

# Create color vector for row names based on pathway classification
row_colors <- sapply(row_names, get_pathway_color)

# Create black color vector for column names
col_colors <- rep("black", length(col_names))


# Create color vector for row names based on pathway classification
row_colors <- sapply(row_names, get_pathway_color)

pdf(paste0('/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/', cell_text, '_cell_modulecomparisons_allparticipants_wPET_spearman.pdf'), 
    width = 20, height = 20)

corrplot(corr_subset, 
         method = "color",
         p.mat = p_subset,           
         sig.level = 0.05,           
         insig = "label_sig",        
         number.cex = 1.2,           
         tl.cex = 1.5,
         tl.col = "black",           # Set all labels to black first
         cl.cex = 1.2)              

dev.off()

png(paste0('/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/', cell_text, '_modulecomparisons_allparticipants_wPET_spearman.png'), 
    width = 20, height = 20, units = 'in', res = 300)
corrplot(corr_subset, 
         method = "color",
         p.mat = p_subset,           
         sig.level = 0.05,           
         insig = "label_sig",        
         number.cex = 1.2,           
         tl.cex = 1.5,
         tl.col = "black",           # Set all labels to black first
         cl.cex = 1.2)              

dev.off()






library(corrplot)
library(dplyr)


# Check sample size
cat("Number of T2D participants:", nrow(t2d_only), "\n")

# Select relevant columns (remove ID columns)
t2d_data <- t2d_only %>% 
  dplyr::select(-record_id, -mrn, -group)

# Calculate correlations
t2d_corr <- cor(t2d_data, use = 'pairwise.complete.obs', method = 'spearman')

# Calculate p-values
p_values <- cor.mtest(t2d_data, method = 'spearman')

# Create subset for plotting
corr_subset <- as.matrix(t2d_corr[c(1:23), c(24:29), drop = FALSE])
p_subset <- as.matrix(p_values$p[c(1:23), c(24:29), drop = FALSE])

# Add meaningful names
row_names <- c('TCA', 'OxPhos', 'Response to Insulin', 
               'Insulin Receptor Signaling', 'Negative Regulation of Insulin Signaling', 
               'Positive Regulation of Insulin Signaling', 'IGF Receptor Binding', 
               'OxPhos (Hallmark)', 'Fatty Acid Metabolism', 'Glycolysis', 
               'mTORC1 Signaling', 'Adipogenesis', 'Peroxisome', 
               'Inflammatory Response', 'IFN Gamma Response', 'IFN Alpha Response', 
               'IL6 JAK/STAT3', 'TNFa NF kB1', 'Complement', 'Allograft Rejection', 
               'Hypoxia', 'ROS Pathway', 'PI3K AKT mTOR1')

col_names <- c('Cortical K2', 'Cortical F', 'Cortical K2/F', 
               'Cortical K2 (voxel)', 'Cortical F (voxel)', 'Cortical K2/F (voxel)')

rownames(corr_subset) <- row_names
colnames(corr_subset) <- col_names
rownames(p_subset) <- row_names
colnames(p_subset) <- col_names

# Check how many significant correlations exist
cat("Number of correlations with p < 0.05:", sum(p_subset < 0.05, na.rm = TRUE), "\n")
cat("Minimum p-value:", min(p_subset, na.rm = TRUE), "\n")

# OPTION 1: Plot without significance markers (cleanest for n=6)
pdf(paste0('/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/', cell_text, '_modulecomparisons_T2Dparticipants_wPET_spearman.pdf'), 
    width = 12, height = 16)

tmp_plot1 <- corrplot(corr_subset, 
         method = "color",
         col = colorRampPalette(c("red", "white", "blue"))(200),
         tl.cex = 1.2,
         tl.col = "black",
         cl.cex = 1.0,
         addCoef.col = "black",  # Add correlation coefficients
         number.cex = 0.8,
         title = "T2D Only: Pathway Correlations with PET Metrics",
         mar = c(0,0,2,0))
print(tmp_plot1)

dev.off()



png(paste0('/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/', cell_text, '_modulecomparisons_T2Dparticipants_wPET_spearman.png'), 
    width = 20, height = 20, units = 'in', res = 300)
tmp_plot1 <- corrplot(corr_subset, 
                      method = "color",
                      col = colorRampPalette(c("red", "white", "blue"))(200),
                      tl.cex = 1.2,
                      tl.col = "black",
                      cl.cex = 1.0,
                      addCoef.col = "black",  # Add correlation coefficients
                      number.cex = 0.8,
                      title = "T2D Only: Pathway Correlations with PET Metrics",
                      mar = c(0,0,2,0))
print(tmp_plot1)


dev.off()



########## Only LC 




library(corrplot)
library(dplyr)


# Check sample size
cat("Number of lc participants:", nrow(lc_only), "\n")

# Select relevant columns (remove ID columns)
lc_data <- lc_only %>% 
  dplyr::select(-record_id, -mrn, -group)

# Calculate correlations
lc_corr <- cor(lc_data, use = 'pairwise.complete.obs', method = 'spearman')

# Calculate p-values
p_values <- cor.mtest(lc_data, method = 'spearman')

# Create subset for plotting
corr_subset <- as.matrix(lc_corr[c(1:23), c(24:29), drop = FALSE])
p_subset <- as.matrix(p_values$p[c(1:23), c(24:29), drop = FALSE])

# Add meaningful names
row_names <- c('TCA', 'OxPhos', 'Response to Insulin', 
               'Insulin Receptor Signaling', 'Negative Regulation of Insulin Signaling', 
               'Positive Regulation of Insulin Signaling', 'IGF Receptor Binding', 
               'OxPhos (Hallmark)', 'Fatty Acid Metabolism', 'Glycolysis', 
               'mTORC1 Signaling', 'Adipogenesis', 'Peroxisome', 
               'Inflammatory Response', 'IFN Gamma Response', 'IFN Alpha Response', 
               'IL6 JAK/STAT3', 'TNFa NF kB1', 'Complement', 'Allograft Rejection', 
               'Hypoxia', 'ROS Pathway', 'PI3K AKT mTOR1')

col_names <- c('Cortical K2', 'Cortical F', 'Cortical K2/F', 
               'Cortical K2 (voxel)', 'Cortical F (voxel)', 'Cortical K2/F (voxel)')

rownames(corr_subset) <- row_names
colnames(corr_subset) <- col_names
rownames(p_subset) <- row_names
colnames(p_subset) <- col_names

# Check how many significant correlations exist
cat("Number of correlations with p < 0.05:", sum(p_subset < 0.05, na.rm = TRUE), "\n")
cat("Minimum p-value:", min(p_subset, na.rm = TRUE), "\n")

# OPTION 1: Plot without significance markers (cleanest for n=6)
pdf(paste0('/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/', cell_text, '_modulecomparisons_LCparticipants_wPET_spearman.pdf'), 
    width = 12, height = 16)

tmp_plot1 <- corrplot(corr_subset, 
                      method = "color",
                      col = colorRampPalette(c("red", "white", "blue"))(200),
                      tl.cex = 1.2,
                      tl.col = "black",
                      cl.cex = 1.0,
                      addCoef.col = "black",  # Add correlation coefficients
                      number.cex = 0.8,
                      title = "LC Only: Pathway Correlations with PET Metrics",
                      mar = c(0,0,2,0))
print(tmp_plot1)

dev.off()



png(paste0('/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/', cell_text, '_modulecomparisons_LCparticipants_wPET_spearman.png'), 
    width = 20, height = 20, units = 'in', res = 300)
tmp_plot1 <- corrplot(corr_subset, 
                      method = "color",
                      col = colorRampPalette(c("red", "white", "blue"))(200),
                      tl.cex = 1.2,
                      tl.col = "black",
                      cl.cex = 1.0,
                      addCoef.col = "black",  # Add correlation coefficients
                      number.cex = 0.8,
                      title = "LC Only: Pathway Correlations with PET Metrics",
                      mar = c(0,0,2,0))
print(tmp_plot1)


dev.off()









############## Test Differences in slope or shift 





remove(list=ls())

module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt')

celltype <- 'PT'



if(celltype == 'PT'){
  module_scores <- module_scores %>% 
    filter(celltype2 == 'PT')
  cell_text <- 'PT'
}else{
  cell_text <- 'All'
}
module_scores$celltype2 <- NULL
module_scores$KPMP_celltype <- NULL





module_scores_summary <- module_scores %>%
  group_by(record_id, mrn, group) %>%
  summarize(
    # Calculate mean for all numeric module score columns
    across(c(TCA_score1, OxPhos_score1, 
             # GO terms
             Response_to_Insulin1, Insulin_Receptor_Signaling1, 
             Neg_Reg_Insulin_Signaling1, Pos_Reg_Insulin_Signaling1, 
             IGF_Receptor_Binding1,
             # HALLMARK pathways - Metabolism
             OxPhos_Hallmark1, Fatty_Acid_Metabolism1, Glycolysis1, 
             mTORC1_Signaling1, Adipogenesis1, Peroxisome1,
             # HALLMARK pathways - Immune/Inflammatory
             Inflammatory_Response1, IFN_Gamma_Response1, IFN_Alpha_Response1,
             IL6_JAK_STAT31, TNFa_NF_kB1, Complement1, Allograft_Rejection1,
             # HALLMARK pathways - Cross-cutting
             Hypoxia1, ROS_Pathway1, PI3K_AKT_mTOR1), 
           list(mean = ~mean(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = 'drop'
  )




harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

#date_of_screen
#screen_date

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))

#gbm <- gbm %>% left_join(dat %>% dplyr::select(record_id, mrn), by='record_id')


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

dat_results$avg_c_k2 <- NULL
dat_results$avg_c_f <- NULL

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


# Function to get pathway colors (same as before)
get_pathway_color <- function(pathway_name) {
  custom_scores <- c("TCA", "OxPhos")
  go_terms <- c("Response to Insulin", "Insulin Receptor Signaling", 
                "Negative Regulation of Insulin Signaling", "Positive Regulation of Insulin Signaling", 
                "IGF Receptor Binding")
  metabolism_hallmark <- c("OxPhos (Hallmark)", "Fatty Acid Metabolism", "Glycolysis", 
                           "mTORC1 Signaling", "Adipogenesis", "Peroxisome")
  immune_inflammatory <- c("Inflammatory Response", "IFN Gamma Response", "IFN Alpha Response",
                           "IL6 JAK/STAT3", "TNFa NF kB1", "Complement", "Allograft Rejection")
  cross_cutting <- c("Hypoxia", "ROS Pathway", "PI3K AKT mTOR1")
  
  if(pathway_name %in% custom_scores) {
    return("#2E8B57")  # Sea Green
  } else if(pathway_name %in% go_terms) {
    return("#4169E1")  # Royal Blue
  } else if(pathway_name %in% metabolism_hallmark) {
    return("#FF6347")  # Tomato Red
  } else if(pathway_name %in% immune_inflammatory) {
    return("#8A2BE2")  # Blue Violet
  } else if(pathway_name %in% cross_cutting) {
    return("#FF8C00")  # Dark Orange
  } else {
    return("#000000")  # Black for unknown
  }
}

combined_df <- module_scores_summary %>% left_join(dat2 %>% 
                                                     dplyr::select(mrn, avg_c_k2, avg_c_f, avg_c_k2_f, 
                                                                   avg_c_k2_vw, avg_c_f_vw, avg_c_k2_f_vw), by= 'mrn')






comparison_data <- combined_df %>%
  filter(!is.na(avg_c_k2)) %>%
  select(record_id, group, avg_c_k2)

# Statistical test
wilcox.test(avg_c_k2 ~ group, data = comparison_data)

# Visualization
ggplot(comparison_data, aes(x = group, y = avg_c_k2, fill = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 3) +
  stat_compare_means(method = "wilcox.test") +
  labs(
       y = "Cortical K2") +
  theme_minimal()


comparison_data <- combined_df %>%
  filter(!is.na(avg_c_k2)) %>%
  select(record_id, group, avg_c_k2_f)

# Statistical test
wilcox.test(avg_c_k2_f ~ group, data = comparison_data)

# Visualization
ggplot(comparison_data, aes(x = group, y = avg_c_k2_f, fill = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 3) +
  stat_compare_means(method = "wilcox.test") +
  labs(
       y = "Cortical K2/f") +
  theme_minimal()




















































########### Only T2D 

combined_df <- module_scores_summary %>% left_join(dat2 %>% 
                                                     dplyr::select(record_id, avg_c_k2, avg_c_f, avg_c_k2_f, 
                                                                   avg_c_k2_vw, avg_c_f_vw, avg_c_k2_f_vw)) %>% 
  filter(group == 'Type_2_Diabetes')

combined_df <- combined_df %>% dplyr::select(-record_id, -mrn, -group) %>%
  filter(!is.na(avg_c_k2))

library(corrplot)

# Calculate correlations
combined_df_corr <- cor(combined_df, use = 'pairwise.complete.obs', method = 'spearman')

# Calculate p-values using cor.mtest
p_values <- cor.mtest(combined_df, method = 'spearman')

# Create subset for plotting
corr_subset <- as.matrix(combined_df_corr[c(1:23), c(24:29), drop = F])
p_subset <- as.matrix(p_values$p[c(1:23), c(24:29), drop = F])

# Add meaningful names
row_names <- c('TCA', 'OxPhos', 'Response to Insulin', 
               'Insulin Receptor Signaling', 'Negative Regulation of Insulin Signaling', 
               'Positive Regulation of Insulin Signaling', 'IGF Receptor Binding', 
               'OxPhos (Hallmark)', 'Fatty Acid Metabolism', 'Glycolysis', 'mTORC1 Signaling', 
               'Adipogenesis', 'Peroxisome', 'Inflammatory Response', 'IFN Gamma Response', 
               'IFN Alpha Response', 'IL6 JAK/STAT3', 'TNFa NF kB1', 'Complement', 
               'Allograft Rejection', 'Hypoxia', 'ROS Pathway', 'PI3K AKT mTOR1')

col_names <- c('Cortical K2', 'Cortical F', 'Cortical K2/F', 
               'Cortical K2 (voxel)', 'Cortical F (voxel)', 'Cortical K2/F (voxel)')

rownames(corr_subset) <- row_names
colnames(corr_subset) <- col_names

# Apply same names to p-value matrix
rownames(p_subset) <- rownames(corr_subset)
colnames(p_subset) <- colnames(corr_subset)

pdf('/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/modulecomparisons_T2Dparticipants_wPET_spearman.pdf', 
    width = 20, height = 20)

corrplot(corr_subset, 
         method = "color",
         p.mat = p_subset,           
         sig.level = 0.05,           
         insig = "label_sig",        
         number.cex = 1.2,           
         tl.cex = 1.5,
         tl.col = "black",           # Set all labels to black first
         cl.cex = 1.2)              

dev.off()



























######################################################## Investigating UACR in Obese, T2D, and LC 



remove(list=ls())


#gbm <- readxl::read_xlsx("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Pathology_Reports_Morphometrics_Shared/Morphometrics/CHCO Morphometrics update 10-19-23.xlsx")
#gbm <- gbm %>% 
#  dplyr::select(record_id = ID, gbm_thick_arith = `GBM thickness nm (arithmetic mean)`, 
#                gbm_thick_harm = `GBM thickness nm (harmonic mean)`)

#harmonized_data <- read.csv("C:/Users/netio/Documents/Harmonized_data/harmonized_dataset.csv", na = '')

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

#date_of_screen
#screen_date

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))

#gbm <- gbm %>% left_join(dat %>% dplyr::select(record_id, mrn), by='record_id')


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

dat_results$avg_c_k2 <- NULL
dat_results$avg_c_f <- NULL

dat_results <- dat_results %>% bind_cols(tmp_results, tmp_results_vw)


dat_results <- dat_results %>% filter(!is.na(avg_c_k2))

dat_results <- dat_results %>% filter(group %in% c('Lean Control', 'Obese Control', 'Type 2 Diabetes'))
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
dat2$epic_sglti2_1[which(dat2$group == 'Obese Control')] <- 'No'

dat2 <- dat2 %>% filter(epic_sglti2_1 != 'Yes')


# Fix data types before creating the table
library(gtsummary)
library(gt)
library(dplyr)

dat2$group[which(dat2$record_id == 'RH2-39-O')] <- 'Obese Control'

# Convert variables to proper data types
combined_df <- dat2 %>%
  mutate(
    # Ensure continuous variables are numeric
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    hba1c = as.numeric(hba1c),
    
    # Ensure categorical variables are factors or characters
    sex = as.factor(sex),
    race_ethnicity = as.factor(race_ethnicity),
    study = as.factor(study),
    group = as.factor(group),
    epic_sglti2_1 = as.factor(epic_sglti2_1)
  )



# Now create the table with proper data types
desc_table1_fixed <- combined_df %>%
  select(age, sex, race_ethnicity, bmi, hba1c, study, group, epic_sglti2_1) %>%
  tbl_summary(
    by = group,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      hba1c ~ "continuous",
      sex ~ "categorical",
      race_ethnicity ~ "categorical",
      study ~ "categorical",
      epic_sglti2_1 ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      sex ~ "Sex", 
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      study ~ "Study",
      epic_sglti2_1 ~ "SGLT2 Inhibitor Use"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "t.test"
    # Skip categorical p-values if they cause issues
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Group**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

# Save version with epic
desc_table1_fixed %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/demographics_OC_LC_T2D_with_epic_final.png", 
         vwidth = 1200, vheight = 800)







### All groups



combined_df <- dat2 %>% #filter(group %in% c('Type 2 Diabetes', 'Obese Control')) %>% 
  dplyr::select(record_id, avg_c_k2, avg_c_f, avg_c_k2_f, 
                avg_c_k2_vw, avg_c_f_vw, avg_c_k2_f_vw, 
                acr_u)


colSums(is.na(combined_df))

combined_df <- combined_df %>% 
  dplyr::select(avg_c_k2, avg_c_f, avg_c_k2_f, 
                avg_c_k2_vw, avg_c_f_vw, avg_c_k2_f_vw, 
                acr_u)

library(corrplot)

# Calculate correlations
combined_df_corr <- cor(combined_df, use = 'pairwise.complete.obs', method = 'spearman')

# Calculate p-values using cor.mtest
p_values <- cor.mtest(combined_df, method = 'spearman')

# Create subset for plotting
corr_subset <- as.matrix(combined_df_corr[c(7), c(1:6), drop = F])
p_subset <- as.matrix(p_values$p[c(7), c(1:6), drop = F])

# Add meaningful names
rownames(corr_subset) <- c('Urine Albumin-Creatinine Ratio')
colnames(corr_subset) <- c('Cortical K2', 'Cortical F', 'Cortical K2/F', 
                           'Cortical K2 (voxel)', 'Cortical F (voxel)', 'Cortical K2/F (voxel)')

# Apply same names to p-value matrix
rownames(p_subset) <- rownames(corr_subset)
colnames(p_subset) <- colnames(corr_subset)

pdf('/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/LC_OC_T2D_Correlations.pdf', 
    width = 15, height = 5)

# Plot without p-values first
corrplot(corr_subset, 
         method = "color",
         number.cex = 1.2,           
         tl.cex = 1.5,
         tl.col = 'black',
         cl.cex = 1.2)

# Add significance stars manually
sig_stars <- ifelse(p_subset < 0.001, "***",
                    ifelse(p_subset < 0.01, "**",
                           ifelse(p_subset < 0.05, "*", "")))

# Add the stars (you may need to adjust x,y positions)
for(i in 1:6) {
  if(sig_stars[1,i] != "") {
    text(i, 1, sig_stars[1,i], cex = 2, col = "red")
  }
}

dev.off()




#### Obese Controls and T2D 

combined_df <- dat2 %>% filter(group %in% c('Type 2 Diabetes', 'Obese Control')) %>% 
  dplyr::select(record_id, avg_c_k2, avg_c_f, avg_c_k2_f, 
                avg_c_k2_vw, avg_c_f_vw, avg_c_k2_f_vw, 
                acr_u)


colSums(is.na(combined_df))

combined_df <- combined_df %>% 
  dplyr::select(avg_c_k2, avg_c_f, avg_c_k2_f, 
                avg_c_k2_vw, avg_c_f_vw, avg_c_k2_f_vw, 
                acr_u)

library(corrplot)

# Calculate correlations
combined_df_corr <- cor(combined_df, use = 'pairwise.complete.obs', method = 'spearman')

# Calculate p-values using cor.mtest
p_values <- cor.mtest(combined_df, method = 'spearman')

# Create subset for plotting
corr_subset <- as.matrix(combined_df_corr[c(7), c(1:6), drop = F])
p_subset <- as.matrix(p_values$p[c(7), c(1:6), drop = F])

# Add meaningful names
rownames(corr_subset) <- c('Urine Albumin-Creatinine Ratio')
colnames(corr_subset) <- c('Cortical K2', 'Cortical F', 'Cortical K2/F', 
                           'Cortical K2 (voxel)', 'Cortical F (voxel)', 'Cortical K2/F (voxel)')

# Apply same names to p-value matrix
rownames(p_subset) <- rownames(corr_subset)
colnames(p_subset) <- colnames(corr_subset)

pdf('/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/OC_T2D_Correlations.pdf', 
    width = 15, height = 5)

# Plot without p-values first
corrplot(corr_subset, 
         method = "color",
         number.cex = 1.2,           
         tl.cex = 1.5,
         tl.col = 'black',
         cl.cex = 1.2)

# Add significance stars manually
sig_stars <- ifelse(p_subset < 0.001, "***",
                    ifelse(p_subset < 0.01, "**",
                           ifelse(p_subset < 0.05, "*", "")))

# Add the stars (you may need to adjust x,y positions)
for(i in 1:6) {
  if(sig_stars[1,i] != "") {
    text(i, 1, sig_stars[1,i], cex = 2, col = "red")
  }
}

dev.off()


#### Testing for just T2D 

test_df <- dat2 %>% 
  dplyr::select(record_id, group, avg_c_k2, avg_c_f, avg_c_k2_f, 
                avg_c_k2_vw, avg_c_f_vw, avg_c_k2_f_vw, 
                acr_u) %>% 
  filter(group == 'Type 2 Diabetes') %>% 
  dplyr::select(-record_id, -group)

library(corrplot)

# Calculate correlations
combined_df_corr <- cor(test_df, use = 'pairwise.complete.obs', method = 'spearman')

# Calculate p-values using cor.mtest
p_values <- cor.mtest(test_df, method = 'spearman')

# Create subset for plotting
corr_subset <- combined_df_corr[c(7), c(1:6)]
p_subset <- p_values$p[c(7), c(1:6)]

# Add meaningful names
rownames(corr_subset) <- c('Urine Albumin-Creatinine Ratio')
colnames(corr_subset) <- c('Cortical K2', 'Cortical F', 'Cortical K2/F', 
                           'Cortical K2 (voxel)', 'Cortical F (voxel)', 'Cortical K2/F (voxel)')

# Apply same names to p-value matrix
rownames(p_subset) <- rownames(corr_subset)
colnames(p_subset) <- colnames(corr_subset)

pdf('C:/Users/netio/Downloads/Correlations.pdf', width = 20, height = 20)
corrplot(corr_subset, 
         method = "color",
         p.mat = p_subset,           # Add p-values
         sig.level = 0.05,           # Significance level
         insig = "label_sig",        # Show significance markers (* for p<0.05, ** for p<0.01, etc.)
         number.cex = 1.2,           # size of correlation numbers
         tl.cex = 1.5,
         tl.col = 'black',
         cl.cex = 1.2)              # size of color legend
dev.off()




























###### Module score comparisons by group 



module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/GO_pathways_modulescores.txt')
module_scores <- module_scores %>% filter(celltype2 == 'PT')
module_scores$celltype2 <- NULL
module_scores$KPMP_celltype <- NULL


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

dat_results$avg_c_k2 <- NULL
dat_results$avg_m_k2 <- NULL
dat_results$avg_c_f <- NULL


dat_results <- dat_results %>% bind_cols(tmp_results, tmp_results_vw)


#dat_results <- dat_results %>% filter(!is.na(avg_c_k2))

dat_results <- dat_results %>% filter(group %in% c('Lean Control', 'Type 2 Diabetes'))


dat_results <- dat_results %>% 
  dplyr::select(mrn, group, starts_with('avg_c'), age, sex, bmi, hba1c, study, epic_sglti2_1, race_ethnicity)




module_scores_summary <- module_scores %>%
  group_by(mrn) %>%
  summarize(
    # Calculate mean, median, and SD for all numeric columns
    across(c(TCA_score1, OxPhos_score1, Response_to_Insulin1, 
             Insulin_Receptor_Signaling1, Neg_Reg_Insulin_Signaling1, 
             Pos_Reg_Insulin_Signaling1, IGF_Receptor_Binding1,
             Chloride_Homeostasis1, Potassium_Homeostasis1,
             Sodium_Homeostasis1, Anion_Homeostasis1,
             Pos_Reg_Phosphorus_Metabolism1, Pos_Reg_Phosphate_Metabolism1,
             Protein_Catabolic_Regulation1, Renal_Sodium_Transport1,
             Renal_Sodium_Absorption1), 
           list(mean = ~mean(.x, na.rm = TRUE),
                median = ~median(.x, na.rm = TRUE), 
                sd = ~sd(.x, na.rm = TRUE)), 
           .names = "{.col}_{.fn}"),
    .groups = 'drop'
  )


combined_df <- module_scores_summary %>% 
  left_join(dat_results, by=c('mrn'))


library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)
library(broom)

# Function to perform statistical tests and add significance annotations
add_significance <- function(data, variable, group_col = "group") {
  # Filter out NA values in both variable and group columns
  clean_data <- data %>% filter(!is.na(.data[[group_col]]) & !is.na(.data[[variable]]))
  
  # Check if we have enough data points
  if(nrow(clean_data) < 3) {
    return(list(significant = FALSE, p_value = NA))
  }
  
  # Get unique groups
  groups <- unique(clean_data[[group_col]])
  
  # If only one group, can't do statistical test
  if(length(groups) < 2) {
    return(list(significant = FALSE, p_value = NA))
  }
  
  # If exactly two groups, perform t-test
  if(length(groups) == 2) {
    group1_data <- clean_data[clean_data[[group_col]] == groups[1], variable]
    group2_data <- clean_data[clean_data[[group_col]] == groups[2], variable]
    
    # Perform unpaired t-test (assuming unequal variances)
    t_result <- t.test(group1_data, group2_data, var.equal = FALSE)
    
    return(list(significant = t_result$p.value < 0.05, p_value = t_result$p.value))
  }
  
  # If more than two groups, perform ANOVA (uncorrected)
  if(length(groups) > 2) {
    aov_result <- aov(as.formula(paste(variable, "~", group_col)), data = clean_data)
    p_value <- summary(aov_result)[[1]][["Pr(>F)"]][1]
    
    return(list(significant = p_value < 0.05, p_value = p_value))
  }
}

# Function to create a single boxplot with significance
create_boxplot <- function(data, variable, title = NULL) {
  # Filter out NA values for plotting
  plot_data <- data %>% filter(!is.na(group) & !is.na(.data[[variable]]))
  
  # Check if we have enough data
  if(nrow(plot_data) < 3) {
    return(ggplot() + 
             theme_void() + 
             labs(title = paste("Insufficient data:", title)) +
             theme(plot.title = element_text(size = 10, hjust = 0.5, color = "gray")))
  }
  
  # Get significance results
  sig_results <- add_significance(data, variable)
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = group, y = .data[[variable]], fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 10, hjust = 0.5),
      legend.position = "none"
    ) +
    labs(title = if(is.null(title)) gsub("1_(mean|median|sd)", "", variable) else title,
         y = "Score")+scale_fill_manual(values = c("Lean Control" = "#619CFF", "Type 2 Diabetes" = "#F8766D"))
  
  # Add significance annotation if significant
  if(!is.na(sig_results$p_value) && sig_results$significant) {
    p <- p + annotate("text", x = Inf, y = Inf, 
                      label = paste0("p = ", round(sig_results$p_value, 4), "*"), 
                      hjust = 1.1, vjust = 1.1, size = 3, color = "red")
  }
  
  return(p)
}

# Function to create 4x4 layout for a specific statistic (mean, median, or sd)
create_4x4_layout <- function(data, statistic = "mean") {
  # Define the 16 traits in order
  traits <- c(
    "TCA_score", "OxPhos_score", "Response_to_Insulin", "Insulin_Receptor_Signaling",
    "Neg_Reg_Insulin_Signaling", "Pos_Reg_Insulin_Signaling", "IGF_Receptor_Binding", 
    "Chloride_Homeostasis", "Potassium_Homeostasis", "Sodium_Homeostasis", 
    "Anion_Homeostasis", "Pos_Reg_Phosphorus_Metabolism", "Pos_Reg_Phosphate_Metabolism", 
    "Protein_Catabolic_Regulation", "Renal_Sodium_Transport", "Renal_Sodium_Absorption"
  )
  
  # Create column names for the specific statistic
  columns <- paste0(traits, "1_", statistic)
  
  # Create list of plots - initialize with NULL placeholders
  plots <- vector("list", 16)
  
  for(i in 1:16) {
    if(columns[i] %in% names(data)) {
      plots[[i]] <- create_boxplot(data, columns[i], 
                                   title = gsub("_", " ", gsub("1_.*", "", columns[i])))
    } else {
      # Create an empty placeholder plot for missing variables
      plots[[i]] <- ggplot() + 
        theme_void() + 
        labs(title = paste("Missing:", gsub("_", " ", traits[i]))) +
        theme(plot.title = element_text(size = 10, hjust = 0.5, color = "gray"))
    }
  }
  
  # Remove any remaining NULL elements (safety check)
  plots <- plots[!sapply(plots, is.null)]
  
  # Ensure we have exactly 16 plots
  if(length(plots) < 16) {
    for(i in (length(plots)+1):16) {
      plots[[i]] <- ggplot() + theme_void()
    }
  }
  
  # Combine into 4x4 grid
  combined_plot <- wrap_plots(plots, ncol = 4, nrow = 4)
  
  # Add overall title
  final_plot <- combined_plot + 
    plot_annotation(title = paste("Gene Set Enrichment Scores by Group -", 
                                  toupper(statistic)),
                    theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
  
  return(final_plot)
}

# Create the three 4x4 layouts
# 1. Mean scores
mean_plot <- create_4x4_layout(combined_df, "mean")

# 2. Median scores  
median_plot <- create_4x4_layout(combined_df, "median")

# 3. Standard deviation scores
sd_plot <- create_4x4_layout(combined_df, "sd")

# Display the plots
print(mean_plot)
print(median_plot) 
print(sd_plot)

# Save the plots
ggsave("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/gene_scores_mean_4x4.pdf", mean_plot, width = 16, height = 12, dpi = 300)
ggsave("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/gene_scores_median_4x4.pdf", median_plot, width = 16, height = 12, dpi = 300)
ggsave("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/gene_scores_sd_4x4.pdf", sd_plot, width = 16, height = 12, dpi = 300)



# Save as high-quality PNG files
ggsave("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/gene_scores_mean_4x4.png", mean_plot, width = 16, height = 12, dpi = 300)
ggsave("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/gene_scores_median_4x4.png", median_plot, width = 16, height = 12, dpi = 300)
ggsave("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/gene_scores_sd_4x4.png", sd_plot, width = 16, height = 12, dpi = 300)













##############All participants (same as Ye Ji; UACR)





remove(list=ls())


#gbm <- readxl::read_xlsx("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Pathology_Reports_Morphometrics_Shared/Morphometrics/CHCO Morphometrics update 10-19-23.xlsx")
#gbm <- gbm %>% 
#  dplyr::select(record_id = ID, gbm_thick_arith = `GBM thickness nm (arithmetic mean)`, 
#                gbm_thick_harm = `GBM thickness nm (harmonic mean)`)

#harmonized_data <- read.csv("C:/Users/netio/Documents/Harmonized_data/harmonized_dataset.csv", na = '')

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

#date_of_screen
#screen_date

#dat <- harmonized_data %>% dplyr::select(-dob) %>% 
#  arrange(date_of_screen) %>% 
#  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
#                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
#                   .by = c(record_id, visit))


dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))




#gbm <- gbm %>% left_join(dat %>% dplyr::select(record_id, mrn), by='record_id')


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

dat_results$avg_c_k2 <- NULL
dat_results$avg_c_f <- NULL

dat_results <- dat_results %>% bind_cols(tmp_results, tmp_results_vw)


dat_results <- dat_results %>% filter(!is.na(avg_c_k2))

dat_results <- dat_results %>% filter(group %in% c('Lean Control', 'Obese Control', 'Type 2 Diabetes'))
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
dat2$epic_sglti2_1[which(dat2$group == 'Obese Control')] <- 'No'

#dat2 <- dat2 %>% filter(epic_sglti2_1 != 'Yes')


# Fix data types before creating the table
library(gtsummary)
library(gt)
library(dplyr)

dat2$group[which(dat2$record_id == 'RH2-39-O')] <- 'Obese Control'

dat2 <- dat2 %>% filter(record_id != 'CRC-55')

# Convert variables to proper data types
combined_df <- dat2 %>%
  mutate(
    # Ensure continuous variables are numeric
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    hba1c = as.numeric(hba1c),
    
    # Ensure categorical variables are factors or characters
    sex = as.factor(sex),
    race_ethnicity = as.factor(race_ethnicity),
    study = as.factor(study),
    group = as.factor(group),
    epic_sglti2_1 = as.factor(epic_sglti2_1)
  )



# Now create the table with proper data types
desc_table1_fixed <- combined_df %>%
  select(age, sex, race_ethnicity, bmi, hba1c, study, group, epic_sglti2_1) %>%
  tbl_summary(
    by = group,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      hba1c ~ "continuous",
      sex ~ "categorical",
      race_ethnicity ~ "categorical",
      study ~ "categorical",
      epic_sglti2_1 ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      sex ~ "Sex", 
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      study ~ "Study",
      epic_sglti2_1 ~ "SGLT2 Inhibitor Use"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "t.test"
    # Skip categorical p-values if they cause issues
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Group**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

# Save version with epic
desc_table1_fixed %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/demographics_OC_LC_T2D_with_epic_final.png", 
         vwidth = 1200, vheight = 800)







### All groups



combined_df <- dat2 %>% #filter(group %in% c('Type 2 Diabetes', 'Obese Control')) %>% 
  dplyr::select(record_id, avg_c_k2, avg_c_f, avg_c_k2_f, 
                avg_c_k2_vw, avg_c_f_vw, avg_c_k2_f_vw, 
                acr_u)


colSums(is.na(combined_df))

combined_df <- combined_df %>% 
  dplyr::select(avg_c_k2, avg_c_f, avg_c_k2_f, 
                avg_c_k2_vw, avg_c_f_vw, avg_c_k2_f_vw, 
                acr_u)

library(corrplot)

# Calculate correlations
combined_df_corr <- cor(combined_df, use = 'pairwise.complete.obs', method = 'spearman')

# Calculate p-values using cor.mtest
p_values <- cor.mtest(combined_df, method = 'spearman')

# Create subset for plotting
corr_subset <- as.matrix(combined_df_corr[c(7), c(1:6), drop = F])
p_subset <- as.matrix(p_values$p[c(7), c(1:6), drop = F])


rownames(corr_subset) <- c('Urine Albumin-Creatinine Ratio')
colnames(corr_subset) <- c('Cortical K2', 'Cortical F', 'Cortical K2/F', 
                           'Cortical K2 (voxel)', 'Cortical F (voxel)', 'Cortical K2/F (voxel)')

# Apply same names to p-value matrix
rownames(p_subset) <- rownames(corr_subset)
colnames(p_subset) <- colnames(corr_subset)

pdf('C:/Users/netio/Downloads/OC_LC_T2D_all_Correlations.pdf', width = 20, height = 20)
corrplot(corr_subset, 
         method = "color",
         p.mat = p_subset,           # Add p-values
         sig.level = 0.05,           # Significance level
         insig = "label_sig",        # Show significance markers (* for p<0.05, ** for p<0.01, etc.)
         number.cex = 1.2,           # size of correlation numbers
         tl.cex = 1.5,
         tl.col = 'black',
         cl.cex = 1.2)              # size of color legend
dev.off()











#################################with no calculating of variables 


dat2 <- dat %>% filter(!is.na(avg_c_k2)) %>% 
  filter(group %in% c('Lean Control', 'Obese Control', 'Type 2 Diabetes')) %>% 
  mutate(avg_c_k2_f = avg_c_k2 / avg_c_f)



dat2 <- dat2 %>% filter(record_id != 'CRC-55')
dat2$group[which(dat2$record_id == 'RH2-39-O')] <- 'Obese Control'



combined_df <- dat2 %>% #filter(group %in% c('Type 2 Diabetes', 'Obese Control')) %>% 
  dplyr::select(record_id, avg_c_k2, avg_c_f, avg_c_k2_f, 
#                avg_c_k2_vw, avg_c_f_vw, avg_c_k2_f_vw, 
                acr_u)


colSums(is.na(combined_df))

combined_df <- combined_df %>% 
  dplyr::select(avg_c_k2, avg_c_f, avg_c_k2_f, 
#                avg_c_k2_vw, avg_c_f_vw, avg_c_k2_f_vw, 
                acr_u)

library(corrplot)

# Calculate correlations
combined_df_corr <- cor(combined_df, use = 'pairwise.complete.obs', 
                        method = 'spearman')

# Calculate p-values using cor.mtest
p_values <- cor.mtest(combined_df, method = 'spearman')

# Create subset for plotting
corr_subset <- as.matrix(combined_df_corr[c(4), c(1:3), drop = F])
p_subset <- as.matrix(p_values$p[c(4), c(1:3), drop = F])


rownames(corr_subset) <- c('Urine Albumin-Creatinine Ratio')
colnames(corr_subset) <- c('Cortical K2', 'Cortical F', 'Cortical K2/F')

# Apply same names to p-value matrix
rownames(p_subset) <- rownames(corr_subset)
colnames(p_subset) <- colnames(corr_subset)

pdf('C:/Users/netio/Downloads/OC_LC_T2D_all_use_harmonized_values_Correlations.pdf', width = 20, height = 20)
corrplot(corr_subset, 
         method = "color",
         p.mat = p_subset,           # Add p-values
         sig.level = 0.05,           # Significance level
         insig = "label_sig",        # Show significance markers (* for p<0.05, ** for p<0.01, etc.)
         number.cex = 1.2,           # size of correlation numbers
         tl.cex = 1.5,
         tl.col = 'black',
         cl.cex = 1.2)              # size of color legend
dev.off()











############ Correlate with GBM and arteriosclerosis

library(corrplot)

dat2 <- data.table::fread("/Users/netio/Downloads/UACR_Allparticipants_forGBM.csv")

combined_df <- dat2 %>% #filter(group %in% c('Type 2 Diabetes', 'Obese Control')) %>% 
  dplyr::select(record_id, avg_c_k2, avg_c_f, avg_c_k2_f,
                arteriosclerosis, arteriolohyalinosis, `GBM thickness`, `GBM thickness arith`, `GBM thickness harm`,
                acr_u) %>% 
  mutate(arteriolohyalinosis_severity = ifelse(arteriolohyalinosis == 'no', 0, 
                                               ifelse(arteriolohyalinosis == 'mild', 1, 
                                                      ifelse(arteriolohyalinosis == 'severe', 2, NA))), 
         GBM_thickness = ifelse(`GBM thickness` %in% c('normal', 'normal/thickened?'), 1, 0), 
         arteriosclerosis = ifelse(arteriosclerosis == 'yes', 1,
                                   ifelse(arteriosclerosis == 'no', 0, NA))) %>%
  select(-arteriolohyalinosis, -`GBM thickness`, -record_id)



combined_df_corr <- cor(combined_df, use = 'pairwise.complete.obs', 
                        method = 'spearman')

# Calculate p-values using cor.mtest
p_values <- cor.mtest(combined_df, method = 'spearman')

# Create subset for plotting
corr_subset <- as.matrix(combined_df_corr[c(7, 4,8, 5, 6, 9), c(1:3), drop = F])
p_subset <- as.matrix(p_values$p[c(7, 4,8, 5, 6, 9), c(1:3), drop = F])


rownames(corr_subset) <- c('Urine Albumin-Creatinine Ratio', 'Arteriosclerosis (Y/N)',  'Arteriolohyalinosis Severity', 
                           'GBM Thickness (arith)', 'GBM Thickness (Harm)',
                          'GBM Thickening (Y/N)')
colnames(corr_subset) <- c('Cortical K2', 'Cortical F', 'Cortical K2/F')

# Apply same names to p-value matrix
rownames(p_subset) <- rownames(corr_subset)
colnames(p_subset) <- colnames(corr_subset)

pdf('C:/Users/netio/Downloads/UACR_GBM_correlations.pdf', width = 20, height = 20)
corrplot(corr_subset, 
         method = "color",
         p.mat = p_subset,           # Add p-values
         sig.level = 0.05,           # Significance level
         insig = "label_sig",        # Show significance markers (* for p<0.05, ** for p<0.01, etc.)
         number.cex = 1.2,           # size of correlation numbers
         tl.cex = 1.5,
         tl.col = 'black',
         cl.cex = 1.2)              # size of color legend
dev.off()




#without obese controls. No UACR.  


library(corrplot)

dat2 <- data.table::fread("/Users/netio/Downloads/UACR_Allparticipants_forGBM.csv")

dat2 <- dat2[-str_which(dat2$record_id, pattern = '-O'),]

combined_df <- dat2 %>% #filter(group %in% c('Type 2 Diabetes', 'Obese Control')) %>% 
  dplyr::select(record_id, avg_c_k2, avg_c_f, avg_c_k2_f,
                arteriosclerosis, arteriolohyalinosis, `GBM thickness`, `GBM thickness arith`, `GBM thickness harm`,
                acr_u) %>% 
  mutate(arteriolohyalinosis_severity = ifelse(arteriolohyalinosis == 'no', 0, 
                                               ifelse(arteriolohyalinosis == 'mild', 1, 
                                                      ifelse(arteriolohyalinosis == 'severe', 2, NA))), 
         GBM_thickness = ifelse(`GBM thickness` %in% c('normal', 'normal/thickened?'), 1, 0), 
         arteriosclerosis = ifelse(arteriosclerosis == 'yes', 1,
                                   ifelse(arteriosclerosis == 'no', 0, NA))) %>%
  select(-arteriolohyalinosis, -`GBM thickness`, -record_id)



combined_df_corr <- cor(combined_df, use = 'pairwise.complete.obs', 
                        method = 'spearman')

# Calculate p-values using cor.mtest
p_values <- cor.mtest(combined_df, method = 'spearman')

# Create subset for plotting
corr_subset <- as.matrix(combined_df_corr[c(4,8, 9), c(1:3), drop = F])
p_subset <- as.matrix(p_values$p[c(4,8, 9), c(1:3), drop = F])


rownames(corr_subset) <- c('Arteriosclerosis (Y/N)',  'Arteriolohyalinosis Severity', 
#                           'GBM Thickness (arith)', 'GBM Thickness (Harm)',
                           'GBM Thickening (Y/N)')
colnames(corr_subset) <- c('Cortical K2', 'Cortical F', 'Cortical K2/F')

# Apply same names to p-value matrix
rownames(p_subset) <- rownames(corr_subset)
colnames(p_subset) <- colnames(corr_subset)

pdf('C:/Users/netio/Downloads/UACR_GBM_correlations_noOC.pdf', width = 20, height = 20)
corrplot(corr_subset, 
         method = "color",
         p.mat = p_subset,           # Add p-values
         sig.level = 0.05,           # Significance level
         insig = "label_sig",        # Show significance markers (* for p<0.05, ** for p<0.01, etc.)
         number.cex = 1.2,           # size of correlation numbers
         tl.cex = 1.5,
         tl.col = 'black',
         cl.cex = 1.2)              # size of color legend
dev.off()


png('C:/Users/netio/Downloads/UACR_GBM_correlations_noOC.png', width = 25, height = 20, units = 'in', res = 300)
corrplot(corr_subset, 
         method = "color",
         p.mat = p_subset,           # Add p-values
         sig.level = 0.05,           # Significance level
         insig = "label_sig",        # Show significance markers (* for p<0.05, ** for p<0.01, etc.)
         number.cex = 2,           # size of correlation numbers
         tl.cex = 3,
         tl.col = 'black',
         cl.cex = 3)              # size of color legend
dev.off()













