########Make Publication-level figures for ROCKIES 





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
  gtsave("C:/Users/netio/Documents/UofW/Rockies/publication_figures/demographics_PET_with_epic_final.png", 
         vwidth = 1200, vheight = 800)







#PET global
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggbreak)

# Prepare data
aim1_df <- dat2 %>% 
  dplyr::select(record_id, group, avg_c_f, avg_c_k2, avg_c_k2_f)

# Reshape data from wide to long format
aim1_long <- aim1_df %>%
  pivot_longer(cols = c(avg_c_f, avg_c_k2, avg_c_k2_f),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = factor(metric, 
                         levels = c("avg_c_f", "avg_c_k2", "avg_c_k2_f"),
                         labels = c("Cortical F", "Cortical K2", "Cortical K2/F")))

base_path <- 'C:/Users/netio/Documents/UofW/Rockies/publication_figures/'

write.table(aim1_long, paste0(base_path, 'PET_figure_aim1_long.txt'), row.names=F, quote=F, sep='\t')




######### Can start here 
aim1_long <- data.table::fread(paste0(base_path, "PET_figure_aim1_long.txt"))
# Calculate statistics for annotations
stat_test <- aim1_long %>%
  group_by(metric) %>%
  pairwise_wilcox_test(value ~ group, p.adjust.method = "none") %>%
  add_significance() %>%
  add_xy_position(x = "metric", dodge = 0.8)

# Adjust statistical test positions for the broken axis
stat_test_adjusted <- stat_test %>%
  mutate(y.position = c(2.9, 0.25, 0.28))

# Create publication-quality base plot
p <- ggplot(aim1_long, aes(x = metric, y = value)) +
  geom_boxplot(aes(fill = group), 
               color = "black", 
               alpha = 0.8, 
               outlier.shape = NA,
               linewidth = 0.5,
               position = position_dodge(width = 0.8)) +
  geom_point(aes(fill = group), 
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.6, 
             size = 3,  # Increased from 2
             shape = 21, 
             color = "black",
             stroke = 0.3) +
  scale_fill_manual(values = c("#c2dfe3", "#fcb1a6"),
                    labels = c("Lean Control", "Type 2 Diabetes")) +
  labs(x = "PET Metrics",
       y = "Value",
       fill = "Group",
       title = "Figure 2.",
       subtitle = "Global PET Imaging Metrics by Group") +
  theme_classic(base_size = 16) +  # Increased from 11
  theme(
    # Text elements
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, color = "black"),  # Increased
    axis.text.y = element_text(size = 20, color = "black"),  # Increased
    axis.title = element_text(size = 22, face = "bold", color = "black"),  # Increased
    
    # Figure label and title
    plot.title = element_text(size = 18, face = "bold", hjust = 0, color = "black"),  # Consistent with other figures
    plot.subtitle = element_text(size = 16, hjust = 0, color = "black",  # Consistent with other figures
                                 margin = margin(b = 10)),
    
    # Legend
    legend.position = "bottom",
    legend.text = element_text(size = 18, color = "black"),  # Increased
    legend.title = element_text(size = 18, face = "bold", color = "black"),  # Increased
    legend.key.size = unit(0.8, "cm"),  # Increased from 0.5
    legend.background = element_rect(fill = "white", color = NA),
    
    # Axis lines
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    
    # Panel
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Add statistical annotations
p_with_stats <- p + 
  stat_pvalue_manual(stat_test_adjusted, 
                     label = "p.adj.signif",
                     tip.length = 0.01,
                     hide.ns = TRUE,
                     size = 6)  # Increased from 4

# Create segmented y-axis plot
p_broken <- p_with_stats + 
  scale_y_break(c(0.3, 0.70), scales = 2)

# Display the plot
print(p_broken)

# Save publication-quality figures with consistent dimensions
ggsave(paste0(base_path, 'Figure2_GlobalPET.pdf'),
       plot = p_broken,
       width = 18,  # Increased to match other figures
       height = 20,  # Increased to match other figures
       units = "in",
       device = cairo_pdf,
       dpi = 300)

ggsave(paste0(base_path, 'Figure2_GlobalPET.tiff'),
       plot = p_broken,
       width = 18,
       height = 20,
       units = "in",
       dpi = 300,
       compression = "lzw")

ggsave(paste0(base_path, 'Figure2_GlobalPET.png'),
       plot = p_broken,
       width = 18,
       height = 20,
       units = "in",
       dpi = 300)











##################################################################################################################################
remove(list=ls())

#Step 4: TCA & OxPhos transcript and pathways in PT cells 

library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(patchwork)  # For combining plots
library(cowplot)    # For extracting legend

base_path <- 'C:/Users/netio/Documents/UofW/Rockies/publication_figures/'

# Lists for storing plots
tca_plots <- list()
oxphos_plots <- list()

lc_files <- list.files('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', pattern='csv')
t2d_files <- list.files('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2/', pattern = 'csv')

celltypes <- c('PT', 'PT-S1/S2', 'PT-S3', 'aPT')

for(i in c(1:length(celltypes))){
  
  celltype <- celltypes[i]
  
  # Set parameters based on cell type
  if(celltype == 'PT'){
    text_size <- 12
    legend_position <- 'bottom'
    axis_angle <- 0
    hjust_val <- 0.5
    asterisk_size <- 8
  } else {
    text_size <- 12
    legend_position <- 'none'
    axis_angle <- 45
    hjust_val <- 1
    asterisk_size <- 6
  }
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  tmp_lc <- lc_files[str_which(lc_files, pattern = paste0('cycle_', celltype2, '_cells'))]
  tmp_t2d <- t2d_files[str_which(t2d_files, pattern = paste0('cycle_', celltype2, '_cells'))]
  
  #========================================
  # TCA
  #========================================
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
  
  if(nrow(tca_full_sig) > 0){
    plot_df <- tca_full_sig %>% 
      gather(key = 'metric_condition', value = 'value', -gene, -direction) %>% 
      separate(metric_condition, into = c('metric', 'condition'), sep = '_') %>% 
      spread(key = metric, value = value) %>% 
      mutate(Significance = ifelse(pvalue < 0.05, '*', ''))
    
    plot_df <- plot_df %>% filter(condition =='lc')
    
    tmp_plot <- ggplot(plot_df, aes(x=gene, y=logFC, fill=direction)) +
      geom_bar(stat='identity', position='dodge', color='black', linewidth=0.3) +
      geom_text(aes(label = Significance), 
                position = position_dodge(width = 0.9),
                vjust = -0.5, 
                size = asterisk_size) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
      labs(x = NULL, 
           y = 'Log2 Fold Change',
           subtitle = celltype) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      scale_fill_manual(name = 'Regulation',
                        labels = c('Up' = 'Up-regulated', 'Down' = 'Down-regulated'),
                        values = c('Down' = '#d73027', 'Up' = '#4575b4')) +
      theme_classic(base_size = 11) +
      theme(
        text = element_text(color = "black"),
        axis.text.x = element_text(angle = axis_angle, hjust = hjust_val, 
                                   size = 10, color = "black",
                                   margin = margin(t = 8)), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11, face = "bold", color = "black"),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks = element_line(linewidth = 0.5, color = "black"),
        legend.position = legend_position,
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        plot.subtitle = element_text(size = 11, hjust = 0.5, face = "bold", 
                                     margin = margin(b = 12)),  # Increased margin
        plot.margin = margin(t = 10, r = 10, b = 20, l = 10),  # Increased margins, especially bottom
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    # Store plot in list
    tca_plots[[celltype]] <- tmp_plot
  }
  
  #========================================
  # OxPhos
  #========================================
  tmp_lc_oxphos <- tmp_lc[str_which(tmp_lc, pattern = 'PHOS_')]
  tmp_t2d_oxphos <- tmp_t2d[str_which(tmp_t2d, pattern = 'PHOS_')]
  
  oxphos_lc <- data.table::fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', tmp_lc_oxphos))
  oxphos_t2d <- data.table::fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2/', tmp_t2d_oxphos))
  
  oxphos_lc <- oxphos_lc %>% dplyr::select(gene, logFC_lc = logFC_groupType_2_Diabetes, 
                                           pvalue_lc = any_of(c("p_groupType_2_Diabetes", "pvalue")))
  
  oxphos_t2d <- oxphos_t2d %>% dplyr::select(gene, logFC_t2d = logFC_epic_sglti2_1Yes, 
                                             pvalue_t2d = any_of(c("p_epic_sglti2_1Yes", "pvalue")))
  
  oxphos_full <- oxphos_lc %>% left_join(oxphos_t2d)
  
  oxphos_full_sig <- oxphos_full %>% 
    mutate(direction = ifelse(logFC_lc < 0, 'Down', 
                              ifelse(logFC_lc > 0, 'Up', NA)))
  
  if(nrow(oxphos_full_sig) > 0){
    plot_df <- oxphos_full_sig %>% 
      gather(key = 'metric_condition', value = 'value', -gene, -direction) %>% 
      separate(metric_condition, into = c('metric', 'condition'), sep = '_') %>% 
      spread(key = metric, value = value) %>% 
      mutate(Significance = ifelse(pvalue < 0.05, '*', ''))
    
    plot_df <- plot_df %>% filter(condition =='lc')
    
    tmp_plot <- ggplot(plot_df, aes(x=gene, y=logFC, fill=direction)) +
      geom_bar(stat='identity', position='dodge', color='black', linewidth=0.3) +
      geom_text(aes(label = Significance), 
                position = position_dodge(width = 0.9),
                vjust = -0.5, 
                size = asterisk_size) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
      labs(x = NULL,
           y = 'Log2 Fold Change',
           subtitle = celltype) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      scale_fill_manual(name = 'Regulation',
                        labels = c('Up' = 'Up-regulated', 'Down' = 'Down-regulated'),
                        values = c('Down' = '#d73027', 'Up' = '#4575b4')) +
      theme_classic(base_size = 11) +
      theme(
        text = element_text(color = "black"),
        axis.text.x = element_text(angle = axis_angle, hjust = hjust_val, 
                                   size = 10, color = "black",
                                   margin = margin(t = 8)),  # Increased margin
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11, face = "bold", color = "black"),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks = element_line(linewidth = 0.5, color = "black"),
        legend.position = legend_position,
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        plot.subtitle = element_text(size = 11, hjust = 0.5, face = "bold",
                                     margin = margin(b = 12)),  # Increased margin
        plot.margin = margin(t = 10, r = 10, b = 20, l = 10),  # Increased margins, especially bottom
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    # Store plot in list
    oxphos_plots[[celltype]] <- tmp_plot
  }
  
  print(paste0(celltype2, ' is done.'))
}

#========================================
# Combine BOTH TCA and OxPhos into Figure 2
#========================================

# Extract legend from the first plot (with legend)
legend <- get_legend(tca_plots[['PT']])

# Remove legends from all plots for cleaner arrangement
tca_plots[['PT']] <- tca_plots[['PT']] + theme(legend.position = "none")
tca_plots[['PT-S1/S2']] <- tca_plots[['PT-S1/S2']] + theme(legend.position = "none")
tca_plots[['PT-S3']] <- tca_plots[['PT-S3']] + theme(legend.position = "none")
tca_plots[['aPT']] <- tca_plots[['aPT']] + theme(legend.position = "none")

oxphos_plots[['PT']] <- oxphos_plots[['PT']] + theme(legend.position = "none")
oxphos_plots[['PT-S1/S2']] <- oxphos_plots[['PT-S1/S2']] + theme(legend.position = "none")
oxphos_plots[['PT-S3']] <- oxphos_plots[['PT-S3']] + theme(legend.position = "none")
oxphos_plots[['aPT']] <- oxphos_plots[['aPT']] + theme(legend.position = "none")

# Add panel labels
tca_plots[['PT']] <- tca_plots[['PT']] + 
  labs(tag = "A") + 
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = c(0.01, 0.98))

tca_plots[['PT-S1/S2']] <- tca_plots[['PT-S1/S2']] + 
  labs(tag = "B") + 
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = c(0.02, 0.98))

tca_plots[['PT-S3']] <- tca_plots[['PT-S3']] + 
  labs(tag = "C") + 
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = c(0.02, 0.98))

tca_plots[['aPT']] <- tca_plots[['aPT']] + 
  labs(tag = "D") + 
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = c(0.02, 0.98))

oxphos_plots[['PT']] <- oxphos_plots[['PT']] + 
  labs(tag = "E") + 
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = c(0.01, 0.98))

oxphos_plots[['PT-S1/S2']] <- oxphos_plots[['PT-S1/S2']] + 
  labs(tag = "F") + 
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = c(0.02, 0.98))

oxphos_plots[['PT-S3']] <- oxphos_plots[['PT-S3']] + 
  labs(tag = "G") + 
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = c(0.02, 0.98))

oxphos_plots[['aPT']] <- oxphos_plots[['aPT']] + 
  labs(tag = "H") + 
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = c(0.02, 0.98))

# Combine all plots with adjusted heights
all_combined <- (tca_plots[['PT']]) / 
  (tca_plots[['PT-S1/S2']] | tca_plots[['PT-S3']] | tca_plots[['aPT']]) /
  (oxphos_plots[['PT']]) /
  (oxphos_plots[['PT-S1/S2']] | oxphos_plots[['PT-S3']] | oxphos_plots[['aPT']]) +
  plot_layout(heights = c(1.3, 1.3, 1.3, 1.3)) +  # Increased heights
  plot_annotation(
    title = 'Figure 2',
    subtitle = 'TCA cycle and oxidative phosphorylation gene expression in proximal tubule cells',
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0, color = "black"),  # Consistent with other figures
      plot.subtitle = element_text(size = 16, hjust = 0, color = "black",  # Consistent with other figures
                                   margin = margin(b = 10))
    )
  )

# Add shared legend at the bottom (positioned lower to avoid overlap)
all_final <- all_combined + 
  inset_element(legend, left = 0.43, bottom = -0.015, right = 0.57, top = 0.015, align_to = 'full')

# Display the plot
print(all_final)

# Save combined figure with increased width and height
ggsave(paste0(base_path, 'Figure2_TCA_OxPhos_combined.pdf'),
       plot = all_final,
       width = 16,   # Increased from 14 to 16
       height = 22,  # Increased from 18 to 22
       units = "in",
       device = cairo_pdf,
       dpi = 300)

ggsave(paste0(base_path, 'Figure2_TCA_OxPhos_combined.tiff'),
       plot = all_final,
       width = 16,   # Increased from 14 to 16
       height = 22,  # Increased from 18 to 22
       units = "in",
       dpi = 300,
       compression = "lzw")

ggsave(paste0(base_path, 'Figure2_TCA_OxPhos_combined.png'),
       plot = all_final,
       width = 16,   # Increased from 14 to 16
       height = 22,  # Increased from 18 to 22
       units = "in",
       dpi = 300)

print("Figure 2 (combined TCA and OxPhos) saved successfully!")







library(dplyr)
library(stringr)

# Initialize results dataframe
ttest_results <- data.frame(
  Pathway = character(),
  Cell_Type = character(),
  n_genes = integer(),
  Mean_log2FC = numeric(),
  SD_log2FC = numeric(),
  SE_log2FC = numeric(),
  t_statistic = numeric(),
  df = numeric(),
  p_value = numeric(),
  CI_lower = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

lc_files <- list.files('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', pattern='csv')
t2d_files <- list.files('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2/', pattern = 'csv')

# Define cell types - PT overall first, then subtypes
celltypes <- c('PT', 'PT-S1/S2', 'PT-S3', 'aPT')

cat("\n==========================================\n")
cat("Running t-tests for the following cell types:\n")
cat("==========================================\n")
for(i in 1:length(celltypes)) {
  cat(sprintf("%d. %s\n", i, celltypes[i]))
}
cat("\n")

for(i in c(1:length(celltypes))){
  
  celltype <- celltypes[i]
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  tmp_lc <- lc_files[str_which(lc_files, pattern = paste0('cycle_', celltype2, '_cells'))]
  tmp_t2d <- t2d_files[str_which(t2d_files, pattern = paste0('cycle_', celltype2, '_cells'))]
  
  #========================================
  # TCA T-test
  #========================================
  tmp_lc_tca <- tmp_lc[str_which(tmp_lc, pattern = 'TCA')]
  
  tca_lc <- data.table::fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', tmp_lc_tca))
  
  tca_lc <- tca_lc %>% dplyr::select(gene, logFC_lc = logFC_groupType_2_Diabetes, 
                                     pvalue_lc =  any_of(c("p_groupType_2_Diabetes", "pvalue")))
  
  # Extract log2FC values for t-test
  tca_log2fc <- tca_lc$logFC_lc
  
  # Run one-sample t-test
  tca_ttest <- t.test(tca_log2fc, mu = 0, alternative = 'greater')
  
  # Store results
  ttest_results <- rbind(ttest_results, data.frame(
    Pathway = "TCA Cycle",
    Cell_Type = celltype,
    n_genes = length(tca_log2fc),
    Mean_log2FC = mean(tca_log2fc),
    SD_log2FC = sd(tca_log2fc),
    SE_log2FC = sd(tca_log2fc) / sqrt(length(tca_log2fc)),
    t_statistic = as.numeric(tca_ttest$statistic),
    df = as.numeric(tca_ttest$parameter),
    p_value = tca_ttest$p.value,
    CI_lower = tca_ttest$conf.int[1],
    Significant = ifelse(tca_ttest$p.value < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  ))
  
  cat("\n==========================================\n")
  cat(sprintf("%s - TCA Cycle\n", celltype))
  cat("==========================================\n")
  cat(sprintf("Number of genes: %d\n", length(tca_log2fc)))
  cat(sprintf("Mean log2FC: %.4f\n", mean(tca_log2fc)))
  cat(sprintf("SD: %.4f\n", sd(tca_log2fc)))
  cat(sprintf("SE: %.4f\n", sd(tca_log2fc) / sqrt(length(tca_log2fc))))
  cat(sprintf("t-statistic: %.4f\n", tca_ttest$statistic))
  cat(sprintf("df: %d\n", tca_ttest$parameter))
  cat(sprintf("p-value: %.6f\n", tca_ttest$p.value))
  cat(sprintf("95%% CI: (%.4f, Inf)\n", tca_ttest$conf.int[1]))
  cat(sprintf("Significant (p < 0.05): %s\n", ifelse(tca_ttest$p.value < 0.05, "YES", "NO")))
  
  #========================================
  # OxPhos T-test
  #========================================
  tmp_lc_oxphos <- tmp_lc[str_which(tmp_lc, pattern = 'PHOS_')]
  
  oxphos_lc <- data.table::fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', tmp_lc_oxphos))
  
  oxphos_lc <- oxphos_lc %>% dplyr::select(gene, logFC_lc = logFC_groupType_2_Diabetes, 
                                           pvalue_lc = any_of(c("p_groupType_2_Diabetes", "pvalue")))
  
  # Extract log2FC values for t-test
  oxphos_log2fc <- oxphos_lc$logFC_lc
  
  # Run one-sample t-test
  oxphos_ttest <- t.test(oxphos_log2fc, mu = 0, alternative = 'greater')
  
  # Store results
  ttest_results <- rbind(ttest_results, data.frame(
    Pathway = "OxPhos",
    Cell_Type = celltype,
    n_genes = length(oxphos_log2fc),
    Mean_log2FC = mean(oxphos_log2fc),
    SD_log2FC = sd(oxphos_log2fc),
    SE_log2FC = sd(oxphos_log2fc) / sqrt(length(oxphos_log2fc)),
    t_statistic = as.numeric(oxphos_ttest$statistic),
    df = as.numeric(oxphos_ttest$parameter),
    p_value = oxphos_ttest$p.value,
    CI_lower = oxphos_ttest$conf.int[1],
    Significant = ifelse(oxphos_ttest$p.value < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  ))
  
  cat("\n==========================================\n")
  cat(sprintf("%s - OxPhos\n", celltype))
  cat("==========================================\n")
  cat(sprintf("Number of genes: %d\n", length(oxphos_log2fc)))
  cat(sprintf("Mean log2FC: %.4f\n", mean(oxphos_log2fc)))
  cat(sprintf("SD: %.4f\n", sd(oxphos_log2fc)))
  cat(sprintf("SE: %.4f\n", sd(oxphos_log2fc) / sqrt(length(oxphos_log2fc))))
  cat(sprintf("t-statistic: %.4f\n", oxphos_ttest$statistic))
  cat(sprintf("df: %d\n", oxphos_ttest$parameter))
  cat(sprintf("p-value: %.6f\n", oxphos_ttest$p.value))
  cat(sprintf("95%% CI: (%.4f, Inf)\n", oxphos_ttest$conf.int[1]))
  cat(sprintf("Significant (p < 0.05): %s\n", ifelse(oxphos_ttest$p.value < 0.05, "YES", "NO")))
}

cat("\n\n==========================================\n")
cat("SUMMARY TABLE\n")
cat("==========================================\n\n")

# Organize results: PT overall first, then subtypes
ttest_results$Cell_Type <- factor(ttest_results$Cell_Type, 
                                  levels = c('PT', 'PT-S1/S2', 'PT-S3', 'aPT'))
ttest_results <- ttest_results %>%
  arrange(Pathway, Cell_Type)

# Print summary table
print(ttest_results)

# Save results to CSV
write.csv(ttest_results, 
          'C:/Users/netio/Documents/UofW/Rockies/publication_figures/TCA_OxPhos_ttests_all_PT.csv',
          row.names = FALSE)

cat("\n\nResults saved to: TCA_OxPhos_ttests_all_PT.csv\n")

# Create a nicely formatted summary for copying
cat("\n\n==========================================\n")
cat("FORMATTED SUMMARY FOR MANUSCRIPT\n")
cat("==========================================\n\n")

cat("TCA CYCLE:\n")
cat("----------\n")
tca_results <- ttest_results %>% filter(Pathway == "TCA Cycle")
for(i in 1:nrow(tca_results)) {
  row <- tca_results[i,]
  cat(sprintf("%-12s: Mean log2FC = %6.3f ± %5.3f, t(%2d) = %6.3f, p %s\n",
              row$Cell_Type,
              row$Mean_log2FC,
              row$SE_log2FC,
              row$df,
              row$t_statistic,
              ifelse(row$p_value < 0.001, "< 0.001",
                     ifelse(row$p_value < 0.01, sprintf("= %.3f", row$p_value),
                            sprintf("= %.4f", row$p_value)))))
}

cat("\nOXIDATIVE PHOSPHORYLATION:\n")
cat("--------------------------\n")
oxphos_results <- ttest_results %>% filter(Pathway == "OxPhos")
for(i in 1:nrow(oxphos_results)) {
  row <- oxphos_results[i,]
  cat(sprintf("%-12s: Mean log2FC = %6.3f ± %5.3f, t(%2d) = %6.3f, p %s\n",
              row$Cell_Type,
              row$Mean_log2FC,
              row$SE_log2FC,
              row$df,
              row$t_statistic,
              ifelse(row$p_value < 0.001, "< 0.001",
                     ifelse(row$p_value < 0.01, sprintf("= %.3f", row$p_value),
                            sprintf("= %.4f", row$p_value)))))
}

cat("\n==========================================\n")
cat("SUMMARY STATISTICS\n")
cat("==========================================\n")
cat(sprintf("Total comparisons: %d\n", nrow(ttest_results)))
cat(sprintf("Significant results (p < 0.05): %d (%.1f%%)\n", 
            sum(ttest_results$Significant == "Yes"),
            100 * sum(ttest_results$Significant == "Yes") / nrow(ttest_results)))

cat("\nBy Pathway:\n")
cat(sprintf("  TCA Cycle significant: %d/%d\n", 
            sum(tca_results$Significant == "Yes"), 
            nrow(tca_results)))
cat(sprintf("  OxPhos significant: %d/%d\n", 
            sum(oxphos_results$Significant == "Yes"), 
            nrow(oxphos_results)))





### Analyses using GSEA with custom 
library(dplyr)
library(stringr)

# Initialize results dataframe
meta_analysis_results <- data.frame(
  Pathway = character(),
  Cell_Type = character(),
  n_genes = integer(),
  Pooled_log2FC = numeric(),
  Pooled_SE = numeric(),
  z_statistic = numeric(),
  p_value_twosided = numeric(),
  p_value_greater = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

lc_files <- list.files('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', pattern='csv')
t2d_files <- list.files('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2/', pattern = 'csv')

# Define cell types - PT overall first, then subtypes
celltypes <- c('PT', 'PT-S1/S2', 'PT-S3', 'aPT')

cat("\n==========================================\n")
cat("Running Inverse Variance Weighted Meta-Analysis\n")
cat("==========================================\n")
for(i in 1:length(celltypes)) {
  cat(sprintf("%d. %s\n", i, celltypes[i]))
}
cat("\n")

# Function for inverse variance weighted meta-analysis
meta_analysis_geneset <- function(logFC, SE) {
  
  # Inverse variance weights
  weights <- 1 / SE^2
  
  # Pooled effect size
  pooled_logFC <- sum(logFC * weights) / sum(weights)
  pooled_SE <- sqrt(1 / sum(weights))
  
  # Test statistics
  z_stat <- pooled_logFC / pooled_SE
  p_value_twosided <- 2 * pnorm(-abs(z_stat))
  p_value_greater <- pnorm(z_stat, lower.tail = FALSE)  # one-sided: upregulated
  
  # 95% Confidence interval
  CI_lower <- pooled_logFC - 1.96 * pooled_SE
  CI_upper <- pooled_logFC + 1.96 * pooled_SE
  
  return(list(
    n_genes = length(logFC),
    pooled_logFC = pooled_logFC,
    pooled_SE = pooled_SE,
    z_statistic = z_stat,
    p_value_twosided = p_value_twosided,
    p_value_greater = p_value_greater,
    CI_lower = CI_lower,
    CI_upper = CI_upper
  ))
}

for(i in c(1:length(celltypes))){
  
  celltype <- celltypes[i]
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  tmp_lc <- lc_files[str_which(lc_files, pattern = paste0('cycle_', celltype2, '_cells'))]
  tmp_t2d <- t2d_files[str_which(t2d_files, pattern = paste0('cycle_', celltype2, '_cells'))]
  
  #========================================
  # TCA Meta-Analysis
  #========================================
  tmp_lc_tca <- tmp_lc[str_which(tmp_lc, pattern = 'TCA')]
  
  tca_lc <- data.table::fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', tmp_lc_tca))
  
  tca_lc <- tca_lc %>% dplyr::select(
    gene, 
    logFC_lc = logFC_groupType_2_Diabetes,
    SE_lc = se_groupType_2_Diabetes,
    pvalue_lc = any_of(c("p_groupType_2_Diabetes", "pvalue"))
  )
  
  # Run meta-analysis
  tca_meta <- meta_analysis_geneset(
    logFC = tca_lc$logFC_lc,
    SE = tca_lc$SE_lc
  )
  
  # Store results
  meta_analysis_results <- rbind(meta_analysis_results, data.frame(
    Pathway = "TCA Cycle",
    Cell_Type = celltype,
    n_genes = tca_meta$n_genes,
    Pooled_log2FC = tca_meta$pooled_logFC,
    Pooled_SE = tca_meta$pooled_SE,
    z_statistic = tca_meta$z_statistic,
    p_value_twosided = tca_meta$p_value_twosided,
    p_value_greater = tca_meta$p_value_greater,
    CI_lower = tca_meta$CI_lower,
    CI_upper = tca_meta$CI_upper,
    Significant = ifelse(tca_meta$p_value_greater < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  ))
  
  cat("\n==========================================\n")
  cat(sprintf("%s - TCA Cycle Meta-Analysis\n", celltype))
  cat("==========================================\n")
  cat(sprintf("Number of genes: %d\n", tca_meta$n_genes))
  cat(sprintf("Pooled log2FC: %.4f ± %.4f\n", tca_meta$pooled_logFC, tca_meta$pooled_SE))
  cat(sprintf("95%% CI: (%.4f, %.4f)\n", tca_meta$CI_lower, tca_meta$CI_upper))
  cat(sprintf("Z-statistic: %.4f\n", tca_meta$z_statistic))
  cat(sprintf("P-value (two-sided): %.6f\n", tca_meta$p_value_twosided))
  cat(sprintf("P-value (upregulation): %.6f\n", tca_meta$p_value_greater))
  cat(sprintf("Significant upregulation (p < 0.05): %s\n", 
              ifelse(tca_meta$p_value_greater < 0.05, "YES", "NO")))
  
  #========================================
  # OxPhos Meta-Analysis
  #========================================
  tmp_lc_oxphos <- tmp_lc[str_which(tmp_lc, pattern = 'PHOS_')]
  
  oxphos_lc <- data.table::fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', tmp_lc_oxphos))
  
  oxphos_lc <- oxphos_lc %>% dplyr::select(
    gene, 
    logFC_lc = logFC_groupType_2_Diabetes,
    SE_lc = se_groupType_2_Diabetes,
    pvalue_lc = any_of(c("p_groupType_2_Diabetes", "pvalue"))
  )
  
  # Run meta-analysis
  oxphos_meta <- meta_analysis_geneset(
    logFC = oxphos_lc$logFC_lc,
    SE = oxphos_lc$SE_lc
  )
  
  # Store results
  meta_analysis_results <- rbind(meta_analysis_results, data.frame(
    Pathway = "OxPhos",
    Cell_Type = celltype,
    n_genes = oxphos_meta$n_genes,
    Pooled_log2FC = oxphos_meta$pooled_logFC,
    Pooled_SE = oxphos_meta$pooled_SE,
    z_statistic = oxphos_meta$z_statistic,
    p_value_twosided = oxphos_meta$p_value_twosided,
    p_value_greater = oxphos_meta$p_value_greater,
    CI_lower = oxphos_meta$CI_lower,
    CI_upper = oxphos_meta$CI_upper,
    Significant = ifelse(oxphos_meta$p_value_greater < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  ))
  
  cat("\n==========================================\n")
  cat(sprintf("%s - OxPhos Meta-Analysis\n", celltype))
  cat("==========================================\n")
  cat(sprintf("Number of genes: %d\n", oxphos_meta$n_genes))
  cat(sprintf("Pooled log2FC: %.4f ± %.4f\n", oxphos_meta$pooled_logFC, oxphos_meta$pooled_SE))
  cat(sprintf("95%% CI: (%.4f, %.4f)\n", oxphos_meta$CI_lower, oxphos_meta$CI_upper))
  cat(sprintf("Z-statistic: %.4f\n", oxphos_meta$z_statistic))
  cat(sprintf("P-value (two-sided): %.6f\n", oxphos_meta$p_value_twosided))
  cat(sprintf("P-value (upregulation): %.6f\n", oxphos_meta$p_value_greater))
  cat(sprintf("Significant upregulation (p < 0.05): %s\n", 
              ifelse(oxphos_meta$p_value_greater < 0.05, "YES", "NO")))
}

cat("\n\n==========================================\n")
cat("SUMMARY TABLE\n")
cat("==========================================\n\n")

# Organize results: PT overall first, then subtypes
meta_analysis_results$Cell_Type <- factor(meta_analysis_results$Cell_Type, 
                                          levels = c('PT', 'PT-S1/S2', 'PT-S3', 'aPT'))
meta_analysis_results <- meta_analysis_results %>%
  arrange(Pathway, Cell_Type)

# Print summary table
print(meta_analysis_results)

# Save results to CSV
write.csv(meta_analysis_results, 
          'C:/Users/netio/Documents/UofW/Rockies/publication_figures/TCA_OxPhos_MetaAnalysis_all_PT.csv',
          row.names = FALSE)

cat("\n\nResults saved to: TCA_OxPhos_MetaAnalysis_all_PT.csv\n")

# Create a nicely formatted summary for copying
cat("\n\n==========================================\n")
cat("FORMATTED SUMMARY FOR MANUSCRIPT\n")
cat("==========================================\n\n")

cat("TCA CYCLE:\n")
cat("----------\n")
tca_results <- meta_analysis_results %>% filter(Pathway == "TCA Cycle")
for(i in 1:nrow(tca_results)) {
  row <- tca_results[i,]
  cat(sprintf("%-12s: Pooled log2FC = %6.3f ± %5.3f (95%% CI: %.3f to %.3f), Z = %6.3f, p %s\n",
              row$Cell_Type,
              row$Pooled_log2FC,
              row$Pooled_SE,
              row$CI_lower,
              row$CI_upper,
              row$z_statistic,
              ifelse(row$p_value_greater < 0.001, "< 0.001",
                     ifelse(row$p_value_greater < 0.01, sprintf("= %.3f", row$p_value_greater),
                            sprintf("= %.4f", row$p_value_greater)))))
}

cat("\nOXIDATIVE PHOSPHORYLATION:\n")
cat("--------------------------\n")
oxphos_results <- meta_analysis_results %>% filter(Pathway == "OxPhos")
for(i in 1:nrow(oxphos_results)) {
  row <- oxphos_results[i,]
  cat(sprintf("%-12s: Pooled log2FC = %6.3f ± %5.3f (95%% CI: %.3f to %.3f), Z = %6.3f, p %s\n",
              row$Cell_Type,
              row$Pooled_log2FC,
              row$Pooled_SE,
              row$CI_lower,
              row$CI_upper,
              row$z_statistic,
              ifelse(row$p_value_greater < 0.001, "< 0.001",
                     ifelse(row$p_value_greater < 0.01, sprintf("= %.3f", row$p_value_greater),
                            sprintf("= %.4f", row$p_value_greater)))))
}

cat("\n==========================================\n")
cat("SUMMARY STATISTICS\n")
cat("==========================================\n")
cat(sprintf("Total comparisons: %d\n", nrow(meta_analysis_results)))
cat(sprintf("Significant results (p < 0.05): %d (%.1f%%)\n", 
            sum(meta_analysis_results$Significant == "Yes"),
            100 * sum(meta_analysis_results$Significant == "Yes") / nrow(meta_analysis_results)))

cat("\nBy Pathway:\n")
cat(sprintf("  TCA Cycle significant: %d/%d\n", 
            sum(tca_results$Significant == "Yes"), 
            nrow(tca_results)))
cat(sprintf("  OxPhos significant: %d/%d\n", 
            sum(oxphos_results$Significant == "Yes"), 
            nrow(oxphos_results)))













##### PSeudotime


library(slingshot)
#library(condiments)


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
#library(table1)
library(clusterProfiler)
library('org.Hs.eg.db')



setwd('C:/Users/netio/Documents/UofW/Rockies/publication_figures/')
#Get files 

set.seed(123)
load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')


so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

test <- so_subset@meta.data

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


test$epic_sglti2_1 <- NULL

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

dat <- dat %>% 
  dplyr::select(-avg_c_k2, -avg_c_f)

dat <- dat %>% bind_cols(tmp_results)

dat <- dat %>% dplyr::select(record_id, epic_sglti2_1, avg_c_k2, avg_c_k2_f) %>% 
  filter(!is.na(epic_sglti2_1))



test <- test %>% left_join(dat, by='record_id')

so_subset@meta.data$epic_sglti2_1 <- test$epic_sglti2_1
so_subset@meta.data$avg_c_k2 <- test$avg_c_k2
so_subset@meta.data$avg_c_k2_f <- test$avg_c_k2_f
so_subset <- subset(so_subset, epic_sglti2_1 == 'No')


so_subset$celltype1 <- case_when(grepl("PT-",so_subset$celltype_rpca)~"PT",
                                 grepl("TAL-",so_subset$celltype_rpca)~"TAL",
                                 grepl("EC-",so_subset$celltype_rpca)~"EC",
                                 grepl("POD",so_subset$celltype_rpca)~"POD",
                                 grepl("MAC",so_subset$celltype_rpca)~"MAC",
                                 grepl("MON",so_subset$celltype_rpca)~"MON",
                                 grepl("PC-",so_subset$celltype_rpca)~"PC",
                                 grepl("FIB",so_subset$celltype_rpca)~"FIB_MC_VSMC",
                                 grepl("DTL",so_subset$celltype_rpca)~"DTL",
                                 so_subset$celltype_rpca=="DCT"~"DCT",
                                 so_subset$celltype_rpca=="ATL"~"ATL",
                                 so_subset$celltype_rpca=="B"~"B",
                                 so_subset$celltype_rpca=="T"~"T")
so_subset$celltype1 <- as.character(so_subset$celltype1)

so_subset$KPMP_celltype2 <- as.character(so_subset$KPMP_celltype)
so_subset$celltype2 <- ifelse(so_subset$KPMP_celltype=="aPT" | 
                                so_subset$KPMP_celltype=="PT-S1/S2" | 
                                so_subset$KPMP_celltype == "PT-S3","PT",
                              ifelse(grepl("TAL",so_subset$KPMP_celltype),"TAL",
                                     ifelse(grepl("EC-",so_subset$KPMP_celltype),"EC",so_subset$KPMP_celltype2)))


so_subset$DCT_celltype <- ifelse((so_subset$KPMP_celltype=="DCT" | 
                                    so_subset$KPMP_celltype=="dDCT"), "DCT","Non-DCT")


test <- so_subset@meta.data
test2 <- test %>% 
  dplyr::select(record_id, avg_c_k2, avg_c_k2_f, epic_sglti2_1) %>% 
  filter(!duplicated(record_id)) %>% 
  mutate(
    avg_c_k2_f_binary = ifelse(avg_c_k2_f >= median(avg_c_k2_f, na.rm = TRUE), "high", "low"), 
    avg_c_k2_binary = ifelse(avg_c_k2 >= median(avg_c_k2, na.rm= TRUE), 'high', 'low')) %>% 
  dplyr::select(record_id, avg_c_k2_f_binary, avg_c_k2_binary)

test <- test %>% left_join(test2, by='record_id')

so_subset@meta.data$avg_c_k2_f_binary <- test$avg_c_k2_f_binary
so_subset@meta.data$avg_c_k2_binary <- test$avg_c_k2_binary

so_subset <- subset(so_subset, subset = avg_c_k2_f_binary %in% c('high', 'low'))


#Analysis
#PT Cells
so_subset <- subset(so_subset, subset = celltype2 == 'PT')
so_subset <- RunUMAP(so_subset, dims = 1:30)

sling_res <- slingshot(as.SingleCellExperiment(so_subset), clusterLabels = 'KPMP_celltype', 
                       start.clus = 'PT-S1', end.clus = 'aPT', reducedDim = 'UMAP')

so_subset$pseudotime <- slingPseudotime(sling_res)[,1]



dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/publication_figures/'

library(ggplot2)
library(viridis)
library(patchwork)
library(scales)
library(ggridges)

# Extract UMAP coordinates
umap_coords <- Embeddings(so_subset, reduction = "umap")

# Extract slingshot curves
sling_curves <- slingCurves(sling_res)

# Create data frame for plotting
plot_df <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  pseudotime = so_subset$pseudotime,
  celltype = so_subset$KPMP_celltype,  # Using original cell type labels
  avg_c_k2_f_binary = so_subset$avg_c_k2_f_binary
)

# Remove cells with NA pseudotime
plot_df_clean <- plot_df %>% filter(!is.na(pseudotime))

# Plot A: UMAP with pseudotime and cell type labels
p1 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_viridis(option = "plasma", na.value = "grey80") +
  theme_classic() +
  labs(title = "Slingshot Pseudotime on UMAP",
       color = "Pseudotime") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Add slingshot curves
for(i in seq_along(sling_curves)) {
  curve_coords <- sling_curves[[i]]$s[sling_curves[[i]]$ord, ]
  p1 <- p1 + geom_path(data = data.frame(UMAP_1 = curve_coords[, 1], 
                                         UMAP_2 = curve_coords[, 2]),
                       aes(x = UMAP_1, y = UMAP_2),
                       color = "black", size = 2, inherit.aes = FALSE)
}

# Add cell type labels at centroids
celltype_centroids <- plot_df %>%
  group_by(celltype) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE)
  )

p1 <- p1 + 
  geom_label(data = celltype_centroids, 
             aes(x = UMAP_1, y = UMAP_2, label = celltype),
             color = "black", fill = "white", alpha = 0.8,
             size = 4, fontface = "bold", inherit.aes = FALSE)

# Plot B: Violin plot comparing pseudotime by Cortical K2/F level
p2 <- ggplot(plot_df_clean, aes(x = avg_c_k2_f_binary, y = pseudotime, fill = avg_c_k2_f_binary)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.1, size = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     label.y = max(plot_df_clean$pseudotime) * 0.99,  # Near top
                     label.x = 1.4) +
  theme_classic() +
  labs(title = "Pseudotime Comparison by Cortical K2/F Level",
       x = "Cortical K2/F Level",
       y = "Pseudotime") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Plot C: Density plot by Cortical K2/F level
p3 <- ggplot(plot_df_clean, aes(x = pseudotime, fill = avg_c_k2_f_binary, color = avg_c_k2_f_binary)) +
  geom_density(alpha = 0.3, size = 1) +
  theme_classic() +
  labs(title = "Cell Density Along Pseudotime Trajectory by Cortical K2/F Level",
       x = "Pseudotime",
       y = "Density",
       fill = "Cortical K2/F",
       color = "Cortical K2/F") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")

# Split pseudotime evenly into 3 regions
min_pseudotime <- min(plot_df_clean$pseudotime, na.rm = TRUE)
max_pseudotime <- max(plot_df_clean$pseudotime, na.rm = TRUE)
pseudotime_range <- max_pseudotime - min_pseudotime

# Define the three regions
region_breaks <- seq(min_pseudotime, max_pseudotime, length.out = 4)
region1_start <- region_breaks[1]
region1_end <- region_breaks[2]
region2_start <- region_breaks[2]
region2_end <- region_breaks[3]
region3_start <- region_breaks[3]
region3_end <- region_breaks[4]

cat("\n=== Pseudotime Regions ===\n")
cat("Region 1 (Early):", round(region1_start, 2), "to", round(region1_end, 2), "\n")
cat("Region 2 (Middle):", round(region2_start, 2), "to", round(region2_end, 2), "\n")
cat("Region 3 (Late):", round(region3_start, 2), "to", round(region3_end, 2), "\n")

# Calculate total cell counts per cell type AND Cortical K2/F level
total_counts_per_celltype_k2f <- plot_df_clean %>%
  count(celltype, avg_c_k2_f_binary, name = "total_cells")

# Create data for all regions and calculate percentages BY Cortical K2/F level
all_region_data <- data.frame()

regions <- list(
  list(name = "Early", start = region1_start, end = region1_end, number = 1),
  list(name = "Middle", start = region2_start, end = region2_end, number = 2),
  list(name = "Late", start = region3_start, end = region3_end, number = 3)
)

for(region in regions) {
  # Extract cells in this region
  region_cells <- plot_df_clean %>%
    filter(pseudotime >= region$start & pseudotime <= region$end)
  
  # Calculate counts in this region BY Cortical K2/F level
  region_counts <- region_cells %>%
    count(celltype, avg_c_k2_f_binary, name = "cells_in_region") %>%
    left_join(total_counts_per_celltype_k2f, by = c("celltype", "avg_c_k2_f_binary")) %>%
    mutate(
      percent_of_celltype = (cells_in_region / total_cells) * 100,
      region_number = region$number,
      region_name = region$name,
      region_start = region$start,
      region_end = region$end
    )
  
  # Add any missing cell type/Cortical K2/F combinations with 0%
  all_combinations <- expand.grid(
    celltype = unique(plot_df_clean$celltype),
    avg_c_k2_f_binary = unique(plot_df_clean$avg_c_k2_f_binary),
    stringsAsFactors = FALSE
  )
  
  existing_combinations <- region_counts %>% 
    select(celltype, avg_c_k2_f_binary) %>%
    distinct()
  
  missing_combinations <- all_combinations %>%
    anti_join(existing_combinations, by = c("celltype", "avg_c_k2_f_binary"))
  
  if(nrow(missing_combinations) > 0) {
    missing_data <- missing_combinations %>%
      left_join(total_counts_per_celltype_k2f, by = c("celltype", "avg_c_k2_f_binary")) %>%
      mutate(
        cells_in_region = 0,
        percent_of_celltype = 0,
        region_number = region$number,
        region_name = region$name,
        region_start = region$start,
        region_end = region$end
      )
    region_counts <- bind_rows(region_counts, missing_data)
  }
  
  all_region_data <- bind_rows(all_region_data, region_counts)
}

# Print detailed summary
cat("\n=== Percentage of Each Cell Type (by Cortical K2/F Level) in Each Region ===\n")
for(i in 1:3) {
  region_data <- all_region_data %>% filter(region_number == i)
  cat("\n", unique(region_data$region_name), "Region (Pseudotime", 
      round(unique(region_data$region_start), 2), "to", 
      round(unique(region_data$region_end), 2), "):\n")
  print(region_data %>% 
          select(celltype, avg_c_k2_f_binary, cells_in_region, total_cells, percent_of_celltype) %>%
          arrange(celltype, avg_c_k2_f_binary))
}

# Save the detailed composition data
write.csv(all_region_data, 
          paste0(dir.results, "Pseudotime_Regions_Celltype_Cortical_K2F_Percentage.csv"), 
          row.names = FALSE)

# Create bar plots showing percentage of each cell type in each region, SPLIT BY Cortical K2/F level
celltype_colors <- c("PT-S1/S2" = "#4DAF4A",  # Green
                     "PT-S3" = "#377EB8",      # Blue
                     "aPT" = "#E41A1C")        # Red
k2f_colors <- c("high" = "#E63946", "low" = "#457B9D")
bar_charts <- list()
for(i in 1:3) {
  region_data <- all_region_data %>% filter(region_number == i)
  region_name <- unique(region_data$region_name)
  region_start <- unique(region_data$region_start)
  region_end <- unique(region_data$region_end)
  
  bar_charts[[i]] <- ggplot(region_data, aes(x = celltype, y = percent_of_celltype, fill = avg_c_k2_f_binary)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.3) +
    scale_fill_manual(values = k2f_colors) +
    geom_text(aes(label = paste0(round(percent_of_celltype, 1), "%")),
              position = position_dodge(width = 0.9),
              vjust = -0.5, size = 3, fontface = "bold") +
    theme_classic() +
    labs(title = paste0(region_name, " Region\n(", round(region_start, 1), " - ", round(region_end, 1), ")"),
         x = "Cell Type",
         y = "% of Cell Type in Region",
         fill = "Cortical K2/F") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
          axis.title = element_text(size = 9),
          legend.position = "right",
          legend.text = element_text(size = 8)) +
    ylim(0, max(all_region_data$percent_of_celltype) * 1.15)  # Add space for labels
}
# Plot G - Ridge plot showing pseudotime distribution by cell type
# Order cell types by median pseudotime
celltype_order <- plot_df_clean %>%
  group_by(celltype) %>%
  summarise(median_pseudotime = median(pseudotime, na.rm = TRUE)) %>%
  arrange(median_pseudotime) %>%
  pull(celltype)
plot_df_clean$celltype_ordered <- factor(plot_df_clean$celltype, 
                                         levels = celltype_order)
p4 <- ggplot(plot_df_clean, aes(x = pseudotime, y = celltype_ordered, fill = celltype_ordered)) +
  geom_density_ridges(alpha = 0.7, scale = 2, rel_min_height = 0.01) +
  scale_fill_manual(values = celltype_colors) +
  # Add vertical lines for region boundaries
  geom_vline(xintercept = c(region1_end, region2_end), linetype = "dashed", color = "black", size = 0.5) +
  theme_classic() +
  labs(title = "Pseudotime Distribution by Cell Type",
       x = "Pseudotime",
       y = "Cell Type") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"))
# Add panel labels
p1_labeled <- p1 + 
  labs(tag = "A") +
  theme(plot.tag = element_text(size = 16, face = "bold"))
p2_labeled <- p2 + 
  labs(tag = "B") +
  theme(plot.tag = element_text(size = 16, face = "bold"))
p3_labeled <- p3 + 
  labs(tag = "C") +
  theme(plot.tag = element_text(size = 16, face = "bold"))
# Label bar charts
for(i in 1:length(bar_charts)) {
  bar_charts[[i]] <- bar_charts[[i]] + 
    labs(tag = LETTERS[3 + i]) +
    theme(plot.tag = element_text(size = 16, face = "bold"))
}
# Label ridge plot as G
p4_labeled <- p4 + 
  labs(tag = "G") +
  theme(plot.tag = element_text(size = 16, face = "bold"))
# Combine all plots
combined_plot <- (p1_labeled | p2_labeled) / 
  (p3_labeled) / 
  (bar_charts[[1]] | bar_charts[[2]] | bar_charts[[3]]) /
  (p4_labeled)
combined_plot <- combined_plot + 
  plot_layout(heights = c(1, 0.8, 0.8, 0.6)) +
  plot_annotation(
    title = "Figure 3.",
    subtitle = 'Pseudotime Analysis in PT Cells by High/Low K2/F',
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
                  plot.subtitle = element_text(size = 16, hjust = 0))
  )
print(combined_plot)
ggsave(paste0(dir.results, "Figure3_Pseudotime_Analysis_Cortical_K2F_with_Ridge.pdf"), 
       plot = combined_plot, width = 16, height = 16)
ggsave(paste0(dir.results, "Figure3_Pseudotime_Analysis_Cortical_K2F_with_Ridge.png"), 
       plot = combined_plot, width = 16, height = 16, dpi = 300)











################ Figure 4 Clinical 




##############All participants (same as Ye Ji; UACR)
remove(list=ls())

dir.results <- "C:/Users/netio/Documents/UofW/Rockies/publication_figures/"


# Figure 4: Clinical Characteristics
library(corrplot)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(patchwork)

# Load data
dat2 <- data.table::fread("/Users/netio/Downloads/UACR_Allparticipants_forGBM.csv")
dat2 <- dat2[-str_which(dat2$record_id, '-O')]

combined_df <- dat2 %>% 
  dplyr::select(record_id, avg_c_k2, avg_c_f, avg_c_k2_f,
                arteriosclerosis, arteriolohyalinosis, `GBM thickness`, 
                `GBM thickness arith`, `GBM thickness harm`,
                acr_u) %>% 
  mutate(arteriolohyalinosis_severity = ifelse(arteriolohyalinosis == 'no', 0, 
                                               ifelse(arteriolohyalinosis == 'mild', 1, 
                                                      ifelse(arteriolohyalinosis == 'severe', 2, NA))), 
         GBM_thickness = `GBM thickness`) %>%
  select(-arteriolohyalinosis, -`GBM thickness`, -record_id)

PET_traits <- c('Cortical K2', 'Cortical F', 'Cortical K2/F')
PET_traits_small <- c('avg_c_k2', 'avg_c_f', 'avg_c_k2_f')

# Define colors: dark blue for 'no', dark red for 'yes'
colors <- c("no" = "#00008B", "yes" = "#8B0000")

# Create GBM thickness boxplots (Panels A-C)
gbm_plots <- list()
for(i in 1:length(PET_traits)){
  tmp_df <- combined_df 
  iter <- which(names(tmp_df) == PET_traits_small[i])
  names(tmp_df)[iter] <- 'Variable'
  
  plot_data <- tmp_df %>% filter(GBM_thickness != '')
  
  gbm_plots[[i]] <- ggplot(plot_data, 
                           aes(x = as.character(GBM_thickness), 
                               y = Variable, 
                               fill = as.character(GBM_thickness))) +
    geom_boxplot() +
    scale_fill_manual(values = colors) +
    stat_compare_means(method = "wilcox.test",
                       label = "p.format",
                       size = 6,
                       p.adjust.method = "none") +
    theme_classic(base_size = 20) +
    theme(
      axis.title = element_text(size = 22, face = "bold"),
      axis.text = element_text(size = 20),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 18),
      plot.tag = element_text(size = 16, face = "bold")
    ) +
    labs(x = 'GBM Thickening', 
         y = PET_traits[i], 
         fill = 'GBM Thickening',
         tag = LETTERS[i])
}

# Create Arteriosclerosis boxplots (Panels D-F)
arterio_plots <- list()
for(i in 1:length(PET_traits)){
  tmp_df <- combined_df 
  iter <- which(names(tmp_df) == PET_traits_small[i])
  names(tmp_df)[iter] <- 'Variable'
  
  plot_data <- tmp_df %>% filter(arteriosclerosis != '')
  
  arterio_plots[[i]] <- ggplot(plot_data, 
                               aes(x = as.character(arteriosclerosis), 
                                   y = Variable, 
                                   fill = as.character(arteriosclerosis))) +
    geom_boxplot() +
    scale_fill_manual(values = colors) +
    stat_compare_means(method = "wilcox.test",
                       label = "p.format",
                       size = 6,
                       p.adjust.method = "none") +
    theme_classic(base_size = 20) +
    theme(
      axis.title = element_text(size = 22, face = "bold"),
      axis.text = element_text(size = 20),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 18),
      plot.tag = element_text(size = 16, face = "bold")
    ) +
    labs(x = 'Arteriosclerosis', 
         y = PET_traits[i], 
         fill = 'Arteriosclerosis',
         tag = LETTERS[3 + i])
}

# Create UACR correlation plot (Panel G)
# Prepare correlation data
dat_uacr <- dat %>% filter(!is.na(avg_c_k2)) %>% 
  filter(group %in% c('Lean Control', 'Obese Control', 'Type 2 Diabetes')) %>% 
  mutate(avg_c_k2_f = avg_c_k2 / avg_c_f)

dat_uacr <- dat_uacr %>% filter(record_id != 'CRC-55')
dat_uacr$group[which(dat_uacr$record_id == 'RH2-39-O')] <- 'Obese Control'

combined_df_uacr <- dat_uacr %>%
  dplyr::select(avg_c_k2, avg_c_f, avg_c_k2_f, acr_u)

# Calculate correlations
combined_df_corr <- cor(combined_df_uacr, use = 'pairwise.complete.obs', 
                        method = 'spearman')

# Calculate p-values
p_values <- cor.mtest(combined_df_uacr, method = 'spearman')

# Create subset for plotting
corr_subset <- as.matrix(combined_df_corr[4, 1:3, drop = FALSE])
p_subset <- as.matrix(p_values$p[4, 1:3, drop = FALSE])

rownames(corr_subset) <- c('UACR')
colnames(corr_subset) <- c('Cortical K2', 'Cortical F', 'Cortical K2/F')
rownames(p_subset) <- rownames(corr_subset)
colnames(p_subset) <- colnames(corr_subset)

# Create correlation plot
temp_file <- paste0(dir.results, "temp_corrplot.png")
png(temp_file, width = 10, height = 4, units = "in", res = 300)
corrplot(corr_subset, 
         method = "color",
         p.mat = p_subset,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 2,
         insig = "label_sig",
         number.cex = 1.5,
         tl.cex = 1.5,
         tl.col = 'black',
         cl.cex = 1.2,
         mar = c(0, 0, 2, 0))
dev.off()

# Convert corrplot to grob for patchwork
corr_grob <- grid::rasterGrob(png::readPNG(temp_file), interpolate = TRUE)
corr_plot <- ggplot() + 
  annotation_custom(corr_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void() +
  labs(tag = "G") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

# Combine all plots
combined_plot <- (gbm_plots[[1]] | gbm_plots[[2]] | gbm_plots[[3]]) /
  (arterio_plots[[1]] | arterio_plots[[2]] | arterio_plots[[3]]) /
  (corr_plot)

combined_plot <- combined_plot + 
  plot_layout(heights = c(1, 1, 0.6)) +
  plot_annotation(
    title = "Figure 4.",
    subtitle = 'Clinical Characteristics and PET Imaging Biomarkers',
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
                  plot.subtitle = element_text(size = 16, hjust = 0))
  )

# Save the figure in all three formats with increased height
ggsave(paste0(dir.results, "Figure4_Clinical_Characteristics.pdf"), 
       plot = combined_plot, width = 18, height = 20)
ggsave(paste0(dir.results, "Figure4_Clinical_Characteristics.png"), 
       plot = combined_plot, width = 18, height = 20, dpi = 300)
ggsave(paste0(dir.results, "Figure4_Clinical_Characteristics.tiff"), 
       plot = combined_plot, width = 18, height = 20, dpi = 300, compression = "lzw")

# Clean up temporary file
file.remove(temp_file)

########## Make demographics for GBM/Arteriosclerosis




library(purrr)
library(corrplot)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(stringr)

dat2 <- data.table::fread("/Users/netio/Downloads/UACR_Allparticipants_forGBM.csv")

dat2 <- dat2[-str_which(dat2$record_id, '-O')]




harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))

dat_small <- dat %>% 
  filter(record_id %in% dat2$record_id) %>% 
  dplyr::select(record_id, age, bmi, hba1c, sex, race_ethnicity, study, group)



dat2 <- dat2 %>% left_join(dat_small, by='record_id')

dat2[dat2 == ""] <- NA

dat2 <- dat2 %>% filter(!is.na(arteriosclerosis))

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
    arteriosclerosis = as.factor(arteriosclerosis),
    arteriolohyalinosis = as.factor(arteriolohyalinosis),
    `GBM thickness` = as.factor(`GBM thickness`)
  ) %>%
  # Filter out rows with missing arteriosclerosis or GBM thickness
  filter(!is.na(arteriosclerosis) & !is.na(`GBM thickness`))

# ===== TABLE 1: Stratified by Arteriosclerosis =====
desc_table_arterio <- combined_df %>%
  select(age, sex, race_ethnicity, bmi, hba1c, study, group, 
         arteriolohyalinosis, `GBM thickness`, arteriosclerosis) %>%
  tbl_summary(
    by = arteriosclerosis,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      hba1c ~ "continuous",
      sex ~ "categorical",
      race_ethnicity ~ "categorical",
      study ~ "categorical",
      group ~ "categorical",
      arteriolohyalinosis ~ 'categorical',
      `GBM thickness` ~ 'categorical'
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
      group ~ "Group",
      arteriolohyalinosis ~ 'Arteriolohyalinosis Severity',
      `GBM thickness` ~ 'GBM Thickening'
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "t.test",
    all_categorical() ~ "fisher.test"
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Arteriosclerosis**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

# Save arteriosclerosis table
desc_table_arterio %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/demographics_arteriosclerosis.png", 
         vwidth = 1400, vheight = 900)

# ===== TABLE 2: Stratified by GBM Thickness =====
desc_table_gbm <- combined_df %>%
  select(age, sex, race_ethnicity, bmi, hba1c, study, group, 
         arteriolohyalinosis, arteriosclerosis, `GBM thickness`) %>%
  tbl_summary(
    by = `GBM thickness`,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      hba1c ~ "continuous",
      sex ~ "categorical",
      race_ethnicity ~ "categorical",
      study ~ "categorical",
      group ~ "categorical",
      arteriolohyalinosis ~ 'categorical',
      arteriosclerosis ~ 'categorical'
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
      group ~ "Group",
      arteriolohyalinosis ~ 'Arteriolohyalinosis Severity',
      arteriosclerosis ~ 'Arteriosclerosis'
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "t.test",
    all_categorical() ~ "fisher.test"
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**GBM Thickening**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

# Save GBM thickness table
desc_table_gbm %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/demographics_GBMthickness.png", 
         vwidth = 1400, vheight = 900)











######### Merge them 

remove(list=ls())
library(qpdf)

dir.results <- '/Users/netio/Documents/UofW/Rockies/publication_figures/'

# First, let's check each PDF for blank pages
pdf_files <- c(
  paste0(dir.results, "Aim1_GlobalPET.pdf"),
  paste0(dir.results, "Figure2_TCA_OxPhos_combined.pdf"),
  paste0(dir.results, "Figure3_Pseudotime_Analysis_Cortical_K2F_with_Ridge.pdf"),
  paste0(dir.results, "Figure4_Clinical_Characteristics.pdf")
)

# Check number of pages in each PDF
for(pdf in pdf_files) {
  info <- pdf_info(pdf)
  cat(basename(pdf), "has", info$pages, "pages\n")
}

# Combine PDFs
pdf_combine(
  input = pdf_files,
  output = paste0(dir.results, "Combined_All_Figures.pdf")
)

# If there's a blank first page in the combined PDF, remove it
pdf_subset(
  input = paste0(dir.results, "Combined_All_Figures.pdf"),
  pages = 2:pdf_length(paste0(dir.results, "Combined_All_Figures.pdf")),
  output = paste0(dir.results, "Combined_All_Figures_NoBlank.pdf")
)












##### Module scores again

library(Seurat)
library(ggpubr)

# Define TCA genes
tca_genes <- c("ACO1", "ACO2", "DLD", "IDH1", "IDH3A", "MDH1", 
               "OGDH", "PDHA1", "SDHD", "SUCLA2", "SUCLG1", "SUCLG2")

# Add TCA module score to your Seurat object
so_subset <- AddModuleScore(
  so_subset,
  features = list(tca_genes),
  name = "TCA_module"
)

# Extract the scores and group info
score_data <- FetchData(so_subset, vars = c("TCA_module1", "group"))

# Statistical test
test_result <- wilcox.test(TCA_module1 ~ group, data = score_data)
print(test_result)

# Get summary stats
score_data %>% 
  group_by(group) %>%
  summarise(
    mean = mean(TCA_module1),
    median = median(TCA_module1),
    sd = sd(TCA_module1)
  )

# Visualize
VlnPlot(so_subset, 
        features = "TCA_module1", 
        group.by = "group",
        pt.size = 0.1) +
  stat_compare_means(method = "wilcox.test", label.y.npc = 0.95) +
  labs(title = "TCA Cycle Module Score",
       y = "Module Score",
       x = "Group") +
  theme_classic()





library(Seurat)
library(ggpubr)
library(dplyr)

# Extract OxPhos genes from your actual data files
# (Use the same files you used in your original script)

# Read the OxPhos file for PT cells (or whichever celltype you want)
celltype <- 'PT'
celltype2 <- str_replace_all(celltype, "/", "_")
celltype2 <- str_replace_all(celltype2, "-", "_")

lc_files <- list.files('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', pattern='csv')
tmp_lc <- lc_files[str_which(lc_files, pattern = paste0('cycle_', celltype2, '_cells'))]
tmp_lc_oxphos <- tmp_lc[str_which(tmp_lc, pattern = 'PHOS_')]

oxphos_data <- data.table::fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', tmp_lc_oxphos))

# Extract the gene list
oxphos_genes <- oxphos_data$gene

print(oxphos_genes)  # Check what genes you have

# Add OxPhos module score to your Seurat object
so_subset <- AddModuleScore(
  so_subset,
  features = list(oxphos_genes),
  name = "OxPhos_module"
)

# Extract the scores and group info
score_data <- FetchData(so_subset, vars = c("OxPhos_module1", "group"))

# Statistical test
test_result <- wilcox.test(OxPhos_module1 ~ group, data = score_data)
print(test_result)

# Get summary stats
score_data %>% 
  group_by(group) %>%
  summarise(
    mean = mean(OxPhos_module1),
    median = median(OxPhos_module1),
    sd = sd(OxPhos_module1)
  )

# Visualize
VlnPlot(so_subset, 
        features = "OxPhos_module1", 
        group.by = "group",
        pt.size = 0.1) +
  stat_compare_means(method = "wilcox.test", label.y.npc = 0.95) +
  labs(title = "Oxidative Phosphorylation Module Score",
       y = "Module Score",
       x = "Group") +
  theme_classic()


##### PET variables


library(corrplot)
library(Seurat)

# Extract module scores and PET variables


harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

dat <- dat %>% 
  dplyr::select(mrn, avg_c_k2, avg_c_f) %>% 
  dplyr::mutate(avg_c_k2_f = avg_c_k2/avg_c_f) %>% 
  dplyr::filter(!is.na(mrn)) %>% 
  dplyr::group_by(mrn) %>% 
  dplyr::summarize(avg_c_k2 = mean(avg_c_k2, na.rm=T), 
                   avg_c_f = mean(avg_c_f, na.rm=T), 
                   avg_c_k2_f = mean(avg_c_k2_f, na.rm=T))

test <- so_subset@meta.data %>% left_join(dat, by='mrn')



so_subset@meta.data$avg_c_k2 <- test$avg_c_k2
so_subset@meta.data$avg_c_f <- test$avg_c_f
so_subset@meta.data$avg_c_k2_f <- test$avg_c_k2_f
correlation_data <- FetchData(so_subset, vars = c("TCA_module1", "OxPhos_module1", 
                                                  "avg_c_k2", "avg_c_f", "avg_c_k2_f"))

# Create correlation matrix: rows = module scores, columns = PET variables
module_scores <- correlation_data[, c("TCA_module1", "OxPhos_module1")]
pet_variables <- correlation_data[, c("avg_c_k2", "avg_c_f", "avg_c_k2_f")]

# Calculate correlation matrix (module scores vs PET variables)
cor_matrix <- cor(module_scores, pet_variables, method = "spearman", use = "pairwise.complete.obs")

# Rename rows for cleaner labels
rownames(cor_matrix) <- c("TCA Module", "OxPhos Module")
colnames(cor_matrix) <- c("Cortical K2", "Cortical F", "Cortical K2/F")

print(cor_matrix)

# Visualize as corrplot
corrplot(cor_matrix, 
         method = "color",
         addCoef.col = "black",
         number.cex = 1.2,
         tl.col = "black",
         tl.cex = 1,
         tl.srt = 45,
         col = colorRampPalette(c("#4575b4", "white", "#d73027"))(200),
         cl.lim = c(-1, 1))

# Alternative: Calculate p-values for each correlation
p_values <- matrix(NA, nrow = 2, ncol = 3)
rownames(p_values) <- c("TCA Module", "OxPhos Module")
colnames(p_values) <- c("Cortical K2", "Cortical F", "Cortical K2/F")

for(i in 1:nrow(module_scores)){
  for(j in 1:ncol(pet_variables)){
    test <- cor.test(module_scores[,i], pet_variables[,j], method = "spearman")
    p_values[i, j] <- test$p.value
  }
}

print("P-values:")
print(p_values)











######Consort Diagrams
################################################################################
# ROCKIES Study CONSORT Diagram with Participant Counts
################################################################################

library(DiagrammeR)
library(dplyr)
library(stringr)
library(purrr)

################################################################################
# Fix Demographics - Handle Co-enrollment by MRN Priority
################################################################################

# Read harmonized data
harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

# First, aggregate by visit/record_id as in your original code
dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

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
dat_results <- dat %>% bind_cols(tmp_results)

# Filter to those with PET data and correct groups
dat_results <- dat_results %>% filter(!is.na(avg_c_k2))
dat_results <- dat_results %>% filter(group %in% c('Lean Control', 'Type 2 Diabetes'))

# Get SGLT2i information (using your original logic)
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

dat2$epic_sglti2_1[which(dat2$group == 'Lean Control')] <- 'No'

################################################################################
# Calculate Participant Counts for CONSORT
################################################################################

# 1. Total enrolled (unique by MRN from full harmonized data)
n_total <- harmonized_data %>% 
  filter(!is.na(mrn)) %>%
  distinct(mrn) %>%
  nrow()

# 2. By group (from full harmonized data)
n_by_group <- harmonized_data %>%
  filter(!is.na(group)) %>%
  group_by(mrn) %>%
  slice(1) %>%
  ungroup() %>%
  count(group)

n_lean <- n_by_group %>% filter(group == "Lean Control") %>% pull(n)
n_obese <- n_by_group %>% filter(group == "Obese Control") %>% pull(n)
n_t2d <- n_by_group %>% filter(group == "Type 2 Diabetes") %>% pull(n)

# Handle missing values
if(length(n_lean) == 0) n_lean <- 0
if(length(n_obese) == 0) n_obese <- 0
if(length(n_t2d) == 0) n_t2d <- 0

# 3. SGLT2i status in T2D (from dat2 which has PET data)
n_t2d_sglt2i_yes <- dat2 %>%
  filter(group == "Type 2 Diabetes", epic_sglti2_1 == "Yes") %>%
  distinct(mrn) %>%
  nrow()

n_t2d_sglt2i_no <- dat2 %>%
  filter(group == "Type 2 Diabetes", (epic_sglti2_1 == "No" | is.na(epic_sglti2_1))) %>%
  distinct(mrn) %>%
  nrow()

# 4. PET Analysis participants (from dat2)
n_pet_lean <- dat2 %>%
  filter(group == "Lean Control") %>%
  distinct(mrn) %>%
  nrow()

n_pet_t2d <- dat2 %>%
  filter(group == "Type 2 Diabetes", (epic_sglti2_1 == "No" | is.na(epic_sglti2_1))) %>%
  distinct(mrn) %>%
  nrow()

n_pet_total <- n_pet_lean + n_pet_t2d

# Final analysis pool (Lean + T2D no SGLT2i)
n_analysis_pool <- n_pet_lean + n_pet_t2d

# 5. scRNA-seq participants (if so_subset exists)
if(exists("so_subset")) {
  n_scrna_total <- so_subset@meta.data %>%
    filter(epic_sglti2_1 == "No" | is.na(epic_sglti2_1)) %>%
    distinct(mrn) %>%
    nrow()
  
  n_scrna_pt <- so_subset@meta.data %>%
    filter((epic_sglti2_1 == "No" | is.na(epic_sglti2_1)), celltype2 == "PT") %>%
    distinct(mrn) %>%
    nrow()
} else {
  n_scrna_total <- "N/A"
  n_scrna_pt <- "N/A"
}

# 6. Clinical characteristics participants
# Read the UACR data
if(file.exists("C:/Users/netio/Downloads/UACR_Allparticipants_forGBM.csv")) {
  dat_clinical <- data.table::fread("C:/Users/netio/Downloads/UACR_Allparticipants_forGBM.csv")
  dat_clinical <- dat_clinical[-str_which(dat_clinical$record_id, '-O')]
  
  n_clinical_total <- dat_clinical %>%
    filter(!is.na(arteriosclerosis) | !is.na(`GBM thickness`)) %>%
    distinct(record_id) %>%
    nrow()
  
  n_clinical_gbm <- dat_clinical %>%
    filter(!is.na(`GBM thickness`), `GBM thickness` != "") %>%
    distinct(record_id) %>%
    nrow()
  
  n_clinical_arterio <- dat_clinical %>%
    filter(!is.na(arteriosclerosis), arteriosclerosis != "") %>%
    distinct(record_id) %>%
    nrow()
} else {
  n_clinical_total <- "N/A"
  n_clinical_gbm <- "N/A"
  n_clinical_arterio <- "N/A"
}

################################################################################
# Create CONSORT Diagram
################################################################################

consort_diagram <- grViz(paste0("
digraph CONSORT {
  
  # Graph attributes
  graph [layout = dot, rankdir = TB, fontsize = 12, splines=ortho, nodesep=0.5, ranksep=0.8]
  
  # Node attributes
  node [shape = box, style = filled, fontname = 'Arial', margin=0.2]
  
  # Define nodes
  enrolled [label = 'Participants Enrolled\\n(Unique by MRN)\\nn = ", n_total, "', 
            fillcolor = '#E8F4F8', width = 3.5, height = 0.8]
  
  grouped [label = 'Grouped by Diagnosis', fillcolor = '#D4E9F7', width = 3, height = 0.6]
  
  lean [label = 'Lean Control\\nn = ", n_lean, "', fillcolor = '#C1E1EC', width = 2.2, height = 0.7]
  obese [label = 'Obese Control\\nn = ", n_obese, "', fillcolor = '#FFE5CC', width = 2.2, height = 0.7]
  t2d [label = 'Type 2 Diabetes\\nn = ", n_t2d, "', fillcolor = '#FFD4D4', width = 2.2, height = 0.7]
  
  exclude_sglt2i [label = 'Excluded:\\nT2D on SGLT2i\\nn = ", n_t2d_sglt2i_yes, "', 
                  fillcolor = '#FFB3B3', style = 'filled,dashed', width = 2.2, height = 0.8]
  
  exclude_obese [label = 'Excluded:\\nObese Control\\nn = ", n_obese, "', 
                 fillcolor = '#FFB3B3', style = 'filled,dashed', width = 2.2, height = 0.8]
  
  analysis_pool [label = 'Available for Analysis\\nLean Control + T2D (no SGLT2i)\\nn = ", 
                                n_analysis_pool, "', fillcolor = '#C8E6C9', width = 4, height = 0.8]
  
  # Analysis branches
  pet_analysis [label = 'PET Imaging Analysis\\n(with PET data)\\nn = ", n_pet_total, 
                                "\\n\\nLean: ", n_pet_lean, " | T2D: ", n_pet_t2d, "', 
                fillcolor = '#E1BEE7', width = 3.2, height = 1.2]
  
  scrna_analysis [label = 'scRNA-seq Analysis\\n(PT cells)\\nn = ", n_scrna_pt, "', 
                  fillcolor = '#BBDEFB', width = 3.2, height = 1]
  
  clinical_analysis [label = 'Clinical Characteristics\\n(Histology + UACR)\\nn = ", 
                                n_clinical_total, 
                                "\\n\\nGBM: ", n_clinical_gbm,
                                " | Arterio: ", n_clinical_arterio, "', 
                     fillcolor = '#FFE082', width = 3.2, height = 1.2]
  
  # Outputs
  fig1 [label = 'Figure 1:\\nGlobal PET Metrics', fillcolor = '#F3E5F5', width = 2.5, height = 0.7]
  fig2 [label = 'Figure 2:\\nTCA/OxPhos Expression', fillcolor = '#E3F2FD', width = 2.5, height = 0.7]
  fig3 [label = 'Figure 3:\\nPseudotime Analysis', fillcolor = '#E3F2FD', width = 2.5, height = 0.7]
  fig4 [label = 'Figure 4:\\nClinical Correlations', fillcolor = '#FFF9C4', width = 2.5, height = 0.7]
  
  # Edges
  enrolled -> grouped
  grouped -> lean
  grouped -> obese
  grouped -> t2d
  
  t2d -> exclude_sglt2i [label = '  SGLT2i = Yes  ', fontsize = 10, style=dashed]
  obese -> exclude_obese [style = dashed]
  
  lean -> analysis_pool
  t2d -> analysis_pool [label = '  SGLT2i = No  ', fontsize = 10]
  
  analysis_pool -> pet_analysis [label = '  Has PET data  ', fontsize = 10]
  analysis_pool -> scrna_analysis [label = '  Has scRNA-seq  ', fontsize = 10]
  analysis_pool -> clinical_analysis [label = '  Has histology  ', fontsize = 10]
  
  pet_analysis -> fig1
  scrna_analysis -> fig2
  scrna_analysis -> fig3
  clinical_analysis -> fig4
  
  # Rank same level nodes together
  { rank = same; lean; obese; t2d }
  { rank = same; exclude_sglt2i; exclude_obese }
  { rank = same; pet_analysis; scrna_analysis; clinical_analysis }
  { rank = same; fig1; fig2; fig3; fig4 }
}
"))

# Display the diagram
print(consort_diagram)

# Save as PDF
consort_diagram %>%
  export_svg() %>%
  charToRaw() %>%
  rsvg::rsvg_pdf(paste0(dir.results, "CONSORT_Diagram.pdf"))

# Save as PNG
consort_diagram %>%
  export_svg() %>%
  charToRaw() %>%
  rsvg::rsvg_png(paste0(dir.results, "CONSORT_Diagram.png"), width = 2400, height = 3000)

################################################################################
# Create Summary Table
################################################################################

consort_summary <- data.frame(
  Stage = c(
    "Total Enrolled (unique MRN)",
    "Lean Control",
    "Obese Control", 
    "Type 2 Diabetes",
    "  - T2D on SGLT2i (excluded)",
    "  - T2D not on SGLT2i",
    "Available for Analysis",
    "PET Imaging Analysis",
    "  - Lean Control",
    "  - Type 2 Diabetes",
    "scRNA-seq Analysis (PT cells)",
    "Clinical Characteristics Analysis"
  ),
  N = c(
    n_total,
    n_lean,
    n_obese,
    n_t2d,
    n_t2d_sglt2i_yes,
    n_t2d_sglt2i_no,
    n_analysis_pool,
    n_pet_total,
    n_pet_lean,
    n_pet_t2d,
    n_scrna_pt,
    n_clinical_total
  )
)

print(consort_summary)

# Save summary table
write.csv(consort_summary, 
          paste0(dir.results, "CONSORT_Summary_Table.csv"),
          row.names = FALSE)

################################################################################
# Print summary to console
################################################################################

cat("\n========================================\n")
cat("ROCKIES STUDY - PARTICIPANT FLOW SUMMARY\n")
cat("========================================\n\n")
cat("Total Enrolled (unique MRN):", n_total, "\n")
cat("  Lean Control:", n_lean, "\n")
cat("  Obese Control:", n_obese, "(excluded from most analyses)\n")
cat("  Type 2 Diabetes:", n_t2d, "\n")
cat("    - On SGLT2i (excluded):", n_t2d_sglt2i_yes, "\n")
cat("    - Not on SGLT2i:", n_t2d_sglt2i_no, "\n\n")
cat("Available for Analysis:", n_analysis_pool, "\n\n")
cat("ANALYSIS-SPECIFIC COUNTS:\n")
cat("  PET Imaging:", n_pet_total, "(Lean:", n_pet_lean, "| T2D:", n_pet_t2d, ")\n")
cat("  scRNA-seq (PT cells):", n_scrna_pt, "\n")
cat("  Clinical Characteristics:", n_clinical_total, "\n")
cat("    - GBM thickness:", n_clinical_gbm, "\n")
cat("    - Arteriosclerosis:", n_clinical_arterio, "\n")
cat("========================================\n\n")






















##### Consort Graph 

########################################
# ROCKIES Study CONSORT Diagram
# Parallel analyses from harmonized dataset
########################################

library(dplyr)
library(ggplot2)
library(DiagrammeR)
library(gtsummary)
library(gt)
library(data.table)
library(stringr)

# Set paths
base_path <- 'C:/Users/netio/Documents/UofW/Rockies/publication_figures/'
harmonized_path <- "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv"

# Load harmonized data
harmonized_data <- read.csv(harmonized_path, na = '')

# ==============================================================================
# STEP 1: Initial Harmonized Dataset
# ==============================================================================
cat("\n========== STEP 1: Initial Harmonized Dataset ==========\n")

dat_initial <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen)

n_initial_total <- nrow(dat_initial)
n_initial_unique <- length(unique(dat_initial$record_id))

cat("Total observations (before filtering):", n_initial_total, "\n")
cat("Unique participants (before filtering):", n_initial_unique, "\n")

# Filter to only include RH, RH2, and CROCODILE studies based on record_id prefix
dat_initial <- dat_initial %>%
  filter(str_detect(record_id, "^RH2|^RH[^2]|^CROC|^CRC"))

n_initial_filtered <- nrow(dat_initial)
n_initial_unique_filtered <- length(unique(dat_initial$record_id))

cat("\nAfter filtering to RH, RH2, CROCODILE:\n")
cat("Total observations:", n_initial_filtered, "\n")
cat("Unique participants:", n_initial_unique_filtered, "\n")
cat("Studies: RENAL-HEIR, RENAL-HEIRitage, CROCODILE\n")

# ==============================================================================
# STEP 2: After aggregating by record_id and visit
# ==============================================================================
cat("\n========== STEP 2: After Data Aggregation ==========\n")

dat <- dat_initial %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

# Standardize study names based on record_id prefix
dat <- dat %>%
  mutate(study = case_when(
    str_detect(record_id, "^RH2") ~ "RENAL-HEIRitage",
    str_detect(record_id, "^RH[^2]") ~ "RENAL-HEIR",
    str_detect(record_id, "^CROC|^CRC") ~ "CROCODILE",
    TRUE ~ study  # Keep existing study name if doesn't match patterns
  ))

n_after_aggregation <- nrow(dat)
cat("After aggregation (by record_id + visit):", n_after_aggregation, "\n")

# Calculate PET averages
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

# Calculate non-volume weighted averages
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

# ==============================================================================
# PARALLEL ANALYSIS 1: FIGURE 1 - Global PET Analysis
# ==============================================================================
cat("\n========== FIGURE 1: Global PET Analysis ==========\n")

# Filter for PET data
dat_pet <- dat_results %>% filter(!is.na(avg_c_k2))
n_with_pet_all <- nrow(dat_pet)

# Filter for LC and T2D
dat_pet <- dat_pet %>% filter(group %in% c('Lean Control', 'Type 2 Diabetes'))
n_pet_lc_t2d <- nrow(dat_pet)

# Determine SGLT2i status
dat_pet$group2 <- NA

RH <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/RENALHEIR-SGLT2.csv')
names(RH) <- c('Subject', 'rep_instr', 'rep_inst', 'SGLT2')

RH2 <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/RenalHEIRitage-SGLT2Use.csv')
names(RH2) <- c('Subject', 'event', 'rep_instr', 'rep_inst', 'mrn', 'SGLT2', 'SGLT2_ever')
RH2 <- RH2 %>% filter(!is.na(mrn))

improve <- data.table::fread('C:/Users/netio/Downloads/IMPROVET2D-SGLT2i_DATA_LABELS_2025-08-25_0938.csv')
names(improve)[5] <- 'SGLT2'
names(improve)[1] <- 'record_id'
improve <- improve %>% filter(!is.na(SGLT2)) %>% filter(SGLT2 != '')

need_med_info <- dat_pet %>% filter(is.na(group2))
improve_small <- improve %>% filter(record_id %in% need_med_info$record_id)
RH_small <- RH %>% filter(Subject %in% need_med_info$record_id)
RH2_small <- RH2 %>% filter(mrn %in% need_med_info$mrn)

for(i in c(1:nrow(RH_small))){
  if(nrow(RH_small) == 0) next
  if(RH_small$SGLT2[i] == 'No'){
    dat_pet$group2[which(dat_pet$record_id == RH_small$Subject[i])] <- 'T2D-No SGLTi2'
    dat_pet$epic_sglti2_1[which(dat_pet$record_id == RH_small$Subject[i])] <- 'No'
  } else if(RH_small$SGLT2[i] == 'Yes'){
    dat_pet$group2[which(dat_pet$record_id == RH_small$Subject[i])] <- 'T2D-SGLTi2'
    dat_pet$epic_sglti2_1[which(dat_pet$record_id == RH_small$Subject[i])] <- 'Yes'
  }
}

for(i in c(1:nrow(RH2_small))){
  if(nrow(RH2_small) == 0) next
  if(RH2_small$SGLT2[i] == 'No'){
    dat_pet$group2[which(dat_pet$mrn == RH2_small$mrn[i])] <- 'T2D-No SGLTi2'
    dat_pet$epic_sglti2_1[which(dat_pet$mrn == RH2_small$mrn[i])] <- 'No'
  } else if(RH2_small$SGLT2[i] == 'Yes'){
    dat_pet$group2[which(dat_pet$mrn == RH2_small$mrn[i])] <- 'T2D-SGLTi2'
    dat_pet$epic_sglti2_1[which(dat_pet$mrn == RH2_small$mrn[i])] <- 'Yes'
  }
}

for(i in c(1:nrow(improve_small))){
  if(nrow(improve_small) == 0) next
  if(improve_small$SGLT2[i] == 'No'){
    dat_pet$group2[which(dat_pet$record_id == improve_small$record_id[i])] <- 'T2D-No SGLTi2'
    dat_pet$epic_sglti2_1[which(dat_pet$record_id == improve_small$record_id[i])] <- 'No'
  } else if(improve_small$SGLT2[i] == 'Yes'){
    dat_pet$group2[which(dat_pet$record_id == improve_small$record_id[i])] <- 'T2D-SGLTi2'
    dat_pet$epic_sglti2_1[which(dat_pet$record_id == improve_small$record_id[i])] <- 'Yes'
  }
}

dat_pet$epic_sglti2_1[which(dat_pet$group == 'Lean Control')] <- 'No'

# Exclude SGLT2i users
dat_fig1 <- dat_pet %>% filter(epic_sglti2_1 != 'Yes')
n_fig1 <- nrow(dat_fig1)
n_fig1_lc <- sum(dat_fig1$group == 'Lean Control')
n_fig1_t2d <- sum(dat_fig1$group == 'Type 2 Diabetes')
n_excluded_sglt2i <- n_pet_lc_t2d - n_fig1

cat("Figure 1 - Final sample:", n_fig1, "\n")
cat("  Lean Control:", n_fig1_lc, "\n")
cat("  Type 2 Diabetes:", n_fig1_t2d, "\n")
cat("  Excluded (SGLT2i):", n_excluded_sglt2i, "\n")

# Calculate demographics for Figure 1
# Lean Control
fig1_lc_demo <- dat_fig1 %>% 
  filter(group == 'Lean Control') %>%
  distinct(record_id, .keep_all = TRUE) %>%
  summarise(
    age_mean = round(mean(age, na.rm=TRUE), 1),
    age_sd = round(sd(age, na.rm=TRUE), 1),
    n_female = sum(sex == 'Female', na.rm=TRUE),
    n_male = sum(sex == 'Male', na.rm=TRUE),
    bmi_mean = round(mean(bmi, na.rm=TRUE), 1),
    bmi_sd = round(sd(bmi, na.rm=TRUE), 1),
    hba1c_mean = round(mean(hba1c, na.rm=TRUE), 2),
    hba1c_sd = round(sd(hba1c, na.rm=TRUE), 2)
  )

fig1_lc_studies <- dat_fig1 %>% 
  filter(group == 'Lean Control') %>%
  distinct(record_id, .keep_all = TRUE) %>%
  group_by(study) %>%
  summarise(n = n()) %>%
  mutate(study_text = paste0(study, ": ", n)) %>%
  pull(study_text) %>%
  paste(collapse = ", ")

# Type 2 Diabetes
fig1_t2d_demo <- dat_fig1 %>% 
  filter(group == 'Type 2 Diabetes') %>%
  distinct(record_id, .keep_all = TRUE) %>%
  summarise(
    age_mean = round(mean(age, na.rm=TRUE), 1),
    age_sd = round(sd(age, na.rm=TRUE), 1),
    n_female = sum(sex == 'Female', na.rm=TRUE),
    n_male = sum(sex == 'Male', na.rm=TRUE),
    bmi_mean = round(mean(bmi, na.rm=TRUE), 1),
    bmi_sd = round(sd(bmi, na.rm=TRUE), 1),
    hba1c_mean = round(mean(hba1c, na.rm=TRUE), 2),
    hba1c_sd = round(sd(hba1c, na.rm=TRUE), 2)
  )

fig1_t2d_studies <- dat_fig1 %>% 
  filter(group == 'Type 2 Diabetes') %>%
  distinct(record_id, .keep_all = TRUE) %>%
  group_by(study) %>%
  summarise(n = n()) %>%
  mutate(study_text = paste0(study, ": ", n)) %>%
  pull(study_text) %>%
  paste(collapse = ", ")

# Demographics Table for Figure 1
combined_df_fig1 <- dat_fig1 %>%
  mutate(
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    hba1c = as.numeric(hba1c),
    sex = as.factor(sex),
    race_ethnicity = as.factor(race_ethnicity),
    study = as.factor(study),
    group = as.factor(group)
  )

desc_table_fig1 <- combined_df_fig1 %>%
  dplyr::select(age, sex, race_ethnicity, bmi, hba1c, study, group) %>%
  tbl_summary(
    by = group,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      hba1c ~ "continuous",
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
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      sex ~ "Sex", 
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
  modify_spanning_header(all_stat_cols() ~ "**Group**")

desc_table_fig1 %>%
  as_gt() %>%
  gtsave(paste0(base_path, "CONSORT_Figure1_Demographics.png"), 
         vwidth = 1200, vheight = 800)

# ==============================================================================
# PARALLEL ANALYSIS 2: FIGURES 2 & 3 - Single-cell RNA-seq
# ==============================================================================
cat("\n========== FIGURES 2 & 3: Single-cell RNA-seq Analysis ==========\n")

# Load single-cell data - CORRECT FILE FOR LC VS T2D
load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')
so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

# Get metadata
sc_meta <- so_subset@meta.data

# Get unique participants in sc-RNA-seq
sc_record_ids <- unique(sc_meta$record_id)
cat("Total unique participants in sc-RNA-seq dataset:", length(sc_record_ids), "\n")

# Get participants from harmonized data
dat_all_participants <- dat_results %>%
  filter(record_id %in% sc_record_ids) %>%
  filter(group %in% c('Lean Control', 'Type 2 Diabetes')) %>%
  dplyr::select(record_id, group, epic_sglti2_1) %>%
  distinct(record_id, .keep_all = TRUE)

cat("\nDEBUG - All sc-RNA-seq participants by group:\n")
print(table(dat_all_participants$group, useNA = "always"))

# Count by group
sc_participants_lc <- dat_all_participants %>%
  filter(group == 'Lean Control')

sc_participants_t2d <- dat_all_participants %>%
  filter(group == 'Type 2 Diabetes' & (epic_sglti2_1 == 'No' | is.na(epic_sglti2_1)))

sc_participants_no_sglt2i <- bind_rows(sc_participants_lc, sc_participants_t2d)

n_sc_lc <- sum(sc_participants_no_sglt2i$group == 'Lean Control', na.rm = TRUE)
n_sc_t2d <- sum(sc_participants_no_sglt2i$group == 'Type 2 Diabetes', na.rm = TRUE)
n_participants_sc <- n_sc_lc + n_sc_t2d

cat("\nParticipants in sc-RNA-seq (no SGLT2i):", n_participants_sc, "\n")
cat("  Lean Control:", n_sc_lc, "\n")
cat("  Type 2 Diabetes:", n_sc_t2d, "\n")

# Figure 3: Subset with BOTH PT cells AND PET data for trajectory analysis
load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')
so_traj <- so_kpmp_sc
remove(so_kpmp_sc)

# Get metadata and filter
traj_meta <- so_traj@meta.data

# Define PT cell types
traj_meta$celltype2 <- ifelse(traj_meta$KPMP_celltype=="aPT" | 
                                traj_meta$KPMP_celltype=="PT-S1/S2" | 
                                traj_meta$KPMP_celltype == "PT-S3", "PT", "other")

# Get PET data
dat_sc_pet <- dat_fig1 %>% 
  dplyr::select(record_id, epic_sglti2_1, avg_c_k2, avg_c_k2_f, group)

# Remove conflicting columns from traj_meta before joining
cols_to_remove_traj <- c('epic_sglti2_1', 'avg_c_k2', 'avg_c_k2_f', 'group')
cols_to_remove_traj <- cols_to_remove_traj[cols_to_remove_traj %in% names(traj_meta)]
if(length(cols_to_remove_traj) > 0) {
  traj_meta <- traj_meta %>% dplyr::select(-all_of(cols_to_remove_traj))
}

# Merge with PET data
traj_with_pet <- traj_meta %>%
  filter(celltype2 == 'PT') %>%
  left_join(dat_sc_pet, by = 'record_id') %>%
  filter(!is.na(avg_c_k2_f)) %>%
  filter(epic_sglti2_1 == 'No' | group == 'Lean Control') %>%
  distinct(record_id) %>%
  pull(record_id)

n_participants_fig3 <- length(traj_with_pet)

cat("\nFigure 3 - Pseudotime (subset with Those with PET data)\n")
cat("  Participants with PT cells AND PET:", n_participants_fig3, "\n")

# Calculate demographics for sc-RNA-seq
# Lean Control
sc_lc_demo <- dat_all_participants %>% 
  filter(group == 'Lean Control') %>%
  distinct(record_id, .keep_all = TRUE) %>%
  left_join(dat_results %>% dplyr::select(record_id, age, sex, bmi, hba1c, study), by = 'record_id') %>%
  summarise(
    age_mean = round(mean(age, na.rm=TRUE), 1),
    age_sd = round(sd(age, na.rm=TRUE), 1),
    n_female = sum(sex == 'Female', na.rm=TRUE),
    n_male = sum(sex == 'Male', na.rm=TRUE),
    bmi_mean = round(mean(bmi, na.rm=TRUE), 1),
    bmi_sd = round(sd(bmi, na.rm=TRUE), 1),
    hba1c_mean = round(mean(hba1c, na.rm=TRUE), 2),
    hba1c_sd = round(sd(hba1c, na.rm=TRUE), 2)
  )

sc_lc_studies <- dat_all_participants %>% 
  filter(group == 'Lean Control') %>%
  distinct(record_id, .keep_all = TRUE) %>%
  left_join(dat_results %>% dplyr::select(record_id, study), by = 'record_id') %>%
  group_by(study) %>%
  summarise(n = n()) %>%
  mutate(study_text = paste0(study, ": ", n)) %>%
  pull(study_text) %>%
  paste(collapse = ", ")

# Type 2 Diabetes (no SGLT2i)
sc_t2d_demo <- sc_participants_t2d %>% 
  distinct(record_id, .keep_all = TRUE) %>%
  left_join(dat_results %>% dplyr::select(record_id, age, sex, bmi, hba1c, study), by = 'record_id') %>%
  summarise(
    age_mean = round(mean(age, na.rm=TRUE), 1),
    age_sd = round(sd(age, na.rm=TRUE), 1),
    n_female = sum(sex == 'Female', na.rm=TRUE),
    n_male = sum(sex == 'Male', na.rm=TRUE),
    bmi_mean = round(mean(bmi, na.rm=TRUE), 1),
    bmi_sd = round(sd(bmi, na.rm=TRUE), 1),
    hba1c_mean = round(mean(hba1c, na.rm=TRUE), 2),
    hba1c_sd = round(sd(hba1c, na.rm=TRUE), 2)
  )

sc_t2d_studies <- sc_participants_t2d %>% 
  distinct(record_id, .keep_all = TRUE) %>%
  left_join(dat_results %>% dplyr::select(record_id, study), by = 'record_id') %>%
  group_by(study) %>%
  summarise(n = n()) %>%
  mutate(study_text = paste0(study, ": ", n)) %>%
  pull(study_text) %>%
  paste(collapse = ", ")

# ==============================================================================
# PARALLEL ANALYSIS 3: FIGURE 4A-F - GBM/Arteriosclerosis
# ==============================================================================
cat("\n========== FIGURE 4A-F: GBM & Arteriosclerosis ==========\n")

dat_gbm <- data.table::fread("/Users/netio/Downloads/UACR_Allparticipants_forGBM.csv")
dat_gbm <- dat_gbm[-str_which(dat_gbm$record_id, '-O')]

# Merge with harmonized data
dat_small <- dat_fig1 %>% 
  dplyr::select(record_id, age, bmi, hba1c, sex, race_ethnicity, study, group)

dat_gbm <- dat_gbm %>% left_join(dat_small, by='record_id')
dat_gbm[dat_gbm == ""] <- NA

# Filter for complete data
dat_gbm_complete <- dat_gbm %>% filter(!is.na(arteriosclerosis) & !is.na(`GBM thickness`))
n_fig4_histology <- nrow(dat_gbm_complete)

# Count by group and histology features
n_gbm_no <- sum(dat_gbm_complete$`GBM thickness` == 'no', na.rm = TRUE)
n_gbm_yes <- sum(dat_gbm_complete$`GBM thickness` == 'yes', na.rm = TRUE)
n_arterio_no <- sum(dat_gbm_complete$arteriosclerosis == 'no', na.rm = TRUE)
n_arterio_yes <- sum(dat_gbm_complete$arteriosclerosis == 'yes', na.rm = TRUE)

cat("Figure 4A-F - Histology subset:", n_fig4_histology, "\n")
cat("  GBM thickening - No:", n_gbm_no, "/ Yes:", n_gbm_yes, "\n")
cat("  Arteriosclerosis - No:", n_arterio_no, "/ Yes:", n_arterio_yes, "\n")

# ==============================================================================
# PARALLEL ANALYSIS 4: FIGURE 4G - UACR Correlation
# ==============================================================================
cat("\n========== FIGURE 4G: UACR Correlation ==========\n")

# UACR analysis includes LC, OC, and T2D
dat_uacr <- dat_results %>% 
  filter(!is.na(avg_c_k2)) %>% 
  filter(group %in% c('Lean Control', 'Obese Control', 'Type 2 Diabetes')) %>%
  filter(!is.na(acr_u))

# Remove specific outliers as in original code
dat_uacr <- dat_uacr %>% filter(record_id != 'CRC-55')
dat_uacr$group[which(dat_uacr$record_id == 'RH2-39-O')] <- 'Obese Control'

n_fig4_uacr <- nrow(dat_uacr)
n_uacr_lc <- sum(dat_uacr$group == 'Lean Control')
n_uacr_oc <- sum(dat_uacr$group == 'Obese Control')
n_uacr_t2d <- sum(dat_uacr$group == 'Type 2 Diabetes')

cat("Figure 4G - UACR correlation:", n_fig4_uacr, "\n")
cat("  Lean Control:", n_uacr_lc, "\n")
cat("  Obese Control:", n_uacr_oc, "\n")
cat("  Type 2 Diabetes:", n_uacr_t2d, "\n")

# ==============================================================================
# CALCULATE OVERALL ENROLLMENT WITH STUDY BREAKDOWN FOR CONSORT
# ==============================================================================
cat("\n========== Calculating Overall Enrollment ==========\n")

# Get unique participants by MRN across all studies
overall_by_mrn <- dat_results %>%
  distinct(mrn, .keep_all = TRUE) %>%
  filter(!is.na(mrn))

n_total_unique_mrn <- nrow(overall_by_mrn)

# Count participants by study (based on record_id after aggregation)
study_counts_by_record <- dat_results %>%
  distinct(record_id, .keep_all = TRUE) %>%
  dplyr::count(study, name = "n_records") %>%
  arrange(desc(n_records))

cat("Total unique participants (by MRN):", n_total_unique_mrn, "\n")
cat("\nParticipants by study (by record_id):\n")
print(study_counts_by_record)

# Create study enrollment text for first box
study_breakdown <- study_counts_by_record %>%
  mutate(study_text = paste0(study, ": n=", n_records)) %>%
  pull(study_text) %>%
  paste(collapse = "\\n")

# ==============================================================================
# CREATE CONSORT DIAGRAM
# ==============================================================================
cat("\n========== Creating CONSORT Diagram ==========\n")

# Build diagram in parts
header <- sprintf("
digraph CONSORT {
  
  graph [layout = dot, rankdir = TB, fontname = Helvetica, splines = ortho, nodesep = 0.8, ranksep = 0.9]
  
  node [shape = box, style = filled, fontname = Helvetica, fontsize = 10]
  
  # Overall enrollment box with study breakdown
  node1 [label = <
    <b>Participants Enrolled (Unique MRN)</b><br/>
    <font point-size='14'><b>n = %s</b></font><br/>
    <br/>
    <font point-size='10'>
    %s<br/>
    </font>
    <br/>
    <font point-size='9'><i>
    Note: Some participants co-enrolled<br/>
    across studies (counted once by MRN)
    </i></font>
  >, fillcolor = '#E8E8F5', penwidth = 2, width = 3.5]
  
  node2 [label = 'After Data Aggregation\\n(by record_id + visit)\\nn = %s', 
         fillcolor = '#E8E8F5']
  
  # Branch point
  branch [label = 'Participants with Available Data', fillcolor = '#F5F5F5', shape = box, style = 'filled,rounded']
  
  # Exclusion note
  excl_note [label = 'SGLT2i users excluded\\nfrom LC/T2D comparisons', 
             fillcolor = '#FFE6E6', shape = note, fontsize = 9]
",
                  n_total_unique_mrn, 
                  gsub("\\\\n", "<br/>", study_breakdown),
                  n_after_aggregation)

# Figure 1 label
fig1_label <- sprintf("
  # PET Branch
  pet1 [label = 'PET Data Available\\nn = %s', fillcolor = '#D4E8F5']
  pet2 [label = 'LC + T2D only\\nn = %s', fillcolor = '#D4E8F5']
  
  fig1 [label = <
    <b>FIGURE 1: Global PET Analysis</b><br/>
    <b>n = %s</b><br/>
    <br/>
    <b>Lean Control (n=%s)</b><br/>
    Age: %.1f±%.1f, Sex: %sF/%sM<br/>
    BMI: %.1f±%.1f, HbA1c: %.2f±%.2f<br/>
    Studies: %s<br/>
    <br/>
    <b>Type 2 Diabetes (n=%s)</b><br/>
    Age: %.1f±%.1f, Sex: %sF/%sM<br/>
    BMI: %.1f±%.1f, HbA1c: %.2f±%.2f<br/>
    Studies: %s
  >, fillcolor = '#C8E6C9', penwidth = 3, width = 3.5]
",
                      n_with_pet_all, n_pet_lc_t2d,
                      n_fig1, 
                      n_fig1_lc, fig1_lc_demo$age_mean, fig1_lc_demo$age_sd, fig1_lc_demo$n_female, fig1_lc_demo$n_male,
                      fig1_lc_demo$bmi_mean, fig1_lc_demo$bmi_sd, fig1_lc_demo$hba1c_mean, fig1_lc_demo$hba1c_sd, 
                      gsub(", ", ", ", fig1_lc_studies),
                      n_fig1_t2d, fig1_t2d_demo$age_mean, fig1_t2d_demo$age_sd, fig1_t2d_demo$n_female, fig1_t2d_demo$n_male,
                      fig1_t2d_demo$bmi_mean, fig1_t2d_demo$bmi_sd, fig1_t2d_demo$hba1c_mean, fig1_t2d_demo$hba1c_sd, 
                      gsub(", ", ", ", fig1_t2d_studies))

# Figures 2-3 labels
sc_label <- sprintf("
  # Single-cell Branch
  sc1 [label = <
    <b>Single-cell RNA-seq</b><br/>
    <b>n = %s</b><br/>
    <br/>
    <b>Lean Control (n=%s)</b><br/>
    Age: %.1f±%.1f, Sex: %sF/%sM<br/>
    BMI: %.1f±%.1f, HbA1c: %.2f±%.2f<br/>
    Studies: %s<br/>
    <br/>
    <b>Type 2 Diabetes (n=%s)</b><br/>
    Age: %.1f±%.1f, Sex: %sF/%sM<br/>
    BMI: %.1f±%.1f, HbA1c: %.2f±%.2f<br/>
    Studies: %s
  >, fillcolor = '#FFF9C4', width = 3.5]
  
  fig2 [label = <
    <b>FIGURE 2</b><br/>
    TCA and OxPhos Expression<br/>
    All cell types<br/>
    <b>n = %s</b>
  >, fillcolor = '#C8E6C9', penwidth = 3]
  
  sc2 [label = 'Subset with PET data\\nn = %s', 
       fillcolor = '#FFF9C4']
       
  fig3 [label = <
    <b>FIGURE 3</b><br/>
    Pseudotime Trajectory<br/>
    PT cells (High vs Low K2/F)<br/>
    <b>n = %s</b>
  >, fillcolor = '#C8E6C9', penwidth = 3]
",
                    n_participants_sc,
                    n_sc_lc, sc_lc_demo$age_mean, sc_lc_demo$age_sd, sc_lc_demo$n_female, sc_lc_demo$n_male,
                    sc_lc_demo$bmi_mean, sc_lc_demo$bmi_sd, sc_lc_demo$hba1c_mean, sc_lc_demo$hba1c_sd, 
                    gsub(", ", ", ", sc_lc_studies),
                    n_sc_t2d, sc_t2d_demo$age_mean, sc_t2d_demo$age_sd, sc_t2d_demo$n_female, sc_t2d_demo$n_male,
                    sc_t2d_demo$bmi_mean, sc_t2d_demo$bmi_sd, sc_t2d_demo$hba1c_mean, sc_t2d_demo$hba1c_sd, 
                    gsub(", ", ", ", sc_t2d_studies),
                    n_participants_sc, n_participants_fig3, n_participants_fig3)

# Figure 4A-F label
fig4ab_label <- sprintf("
  # Histology Branch
  hist1 [label = 'Kidney Biopsy\\nAvailable\\nn = %s', fillcolor = '#FFE8D5']
  
  fig4ab [label = <
    <b>FIGURE 4A-F</b><br/>
    <b>GBM and Arteriosclerosis</b><br/>
    <b>n = %s</b><br/>
    <br/>
    GBM Thickening<br/>
    No: n=%s, Yes: n=%s<br/>
    <br/>
    Arteriosclerosis<br/>
    No: n=%s, Yes: n=%s
  >, fillcolor = '#C8E6C9', penwidth = 3, width = 2.8]
",
                        n_fig4_histology, n_fig4_histology,
                        n_gbm_no, n_gbm_yes, n_arterio_no, n_arterio_yes)

# Figure 4G label
fig4g_label <- sprintf("
  # UACR Branch
  uacr1 [label = 'UACR Data\\nAvailable\\nn = %s', fillcolor = '#E8D5FF']
  
  fig4g [label = <
    <b>FIGURE 4G</b><br/>
    <b>UACR Correlation</b><br/>
    <b>n = %s</b><br/>
    <br/>
    LC: n=%s, OC: n=%s, T2D: n=%s
  >, fillcolor = '#C8E6C9', penwidth = 3, width = 2.8]
",
                       n_fig4_uacr, n_fig4_uacr, n_uacr_lc, n_uacr_oc, n_uacr_t2d)

# Define edges with port positions
# Define edges with improved routing
# Define edges with improved routing
edges <- "
  # Main flow
  node1 -> node2
  node2 -> branch
  
  # Exclusion note placement
  node2 -> excl_note [style=invis]
  
  # PET branch
  branch -> pet1 [label = '  PET  ']
  pet1 -> pet2
  pet2 -> fig1
  
  # Single-cell branch
  branch -> sc1 [label = '  sc-RNA-seq  ']
  sc1 -> fig2
  sc1 -> sc2 [label = '  Subset  ']
  sc2 -> fig3
  
  # Histology branch
  branch -> hist1 [label = '  Biopsy  ']
  hist1 -> fig4ab
  
  # UACR branch
  branch -> uacr1 [label = '  UACR  ']
  uacr1 -> fig4g
  
  # Layout constraints
  {rank = same; pet1; sc1; hist1; uacr1}
  {rank = same; fig1; fig2; fig4ab; fig4g}
}
"

# Combine all parts
consort_text <- paste0(header, fig1_label, sc_label, fig4ab_label, fig4g_label, edges)

# Create and display diagram
final_diagram <- grViz(consort_text)
print(final_diagram)

# Save outputs
library(htmlwidgets)
saveWidget(final_diagram, 
           paste0(base_path, "CONSORT_Diagram.html"), 
           selfcontained = TRUE)

# Try to save as PNG/PDF using webshot if available
if(require(webshot, quietly = TRUE)) {
  webshot(paste0(base_path, "CONSORT_Diagram.html"),
          paste0(base_path, "CONSORT_Diagram.png"),
          vwidth = 1600, vheight = 2000)
  
  webshot(paste0(base_path, "CONSORT_Diagram.html"),
          paste0(base_path, "CONSORT_Diagram.pdf"),
          vwidth = 1600, vheight = 2000)
  
  cat("\n========== CONSORT Diagram Saved ==========\n")
  cat("PNG:", paste0(base_path, "CONSORT_Diagram.png"), "\n")
  cat("PDF:", paste0(base_path, "CONSORT_Diagram.pdf"), "\n")
} else {
  cat("\n========== CONSORT Diagram Saved ==========\n")
  cat("To save as PNG/PDF, install webshot:\n")
  cat("  install.packages('webshot')\n")
  cat("  webshot::install_phantomjs()\n")
}

cat("HTML:", paste0(base_path, "CONSORT_Diagram.html"), "\n")
cat("\n========== CONSORT Diagram Complete ==========\n")






