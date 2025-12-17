## PET global and voxel analysis 
## PET global and voxel analysis 
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

# Modified function to calculate both voxel-weighted and global averages
PET_avg <- function(data){
  # Voxel-weighted calculations
  tmp_df_vw <- data %>% dplyr::select(lc_k2_wo_cyst_vw, rc_k2_wo_cyst_vw, 
                                      lm_k2_wo_cyst_vw, rm_k2_wo_cyst_vw,
                                      lc_f, rc_f, lm_f, rm_f)
  
  avg_c_k2_vw <- tmp_df_vw %>%
    dplyr::select(lc_k2_wo_cyst_vw, rc_k2_wo_cyst_vw) %>% rowMeans(na.rm=T)
  
  avg_m_k2_vw <- tmp_df_vw %>% 
    dplyr::select(lm_k2_wo_cyst_vw, rm_k2_wo_cyst_vw) %>% rowMeans(na.rm=T)
  
  avg_c_f_vw <- tmp_df_vw %>% 
    dplyr::select(lc_f, rc_f) %>% rowMeans(na.rm=T)
  
  avg_m_f_vw <- tmp_df_vw %>% 
    dplyr::select(lm_f, rm_f) %>% rowMeans(na.rm=T)
  
  avg_c_k2_f_vw <- avg_c_k2_vw / avg_c_f_vw
  avg_m_k2_f_vw <- avg_m_k2_vw / avg_m_f_vw
  
  # Global calculations (without _vw suffix)
  tmp_df_global <- data %>% dplyr::select(lc_k2, rc_k2, 
                                          lm_k2, rm_k2,
                                          lc_f, rc_f, 
                                          lm_f, rm_f)
  
  avg_c_k2_global <- tmp_df_global %>%
    dplyr::select(lc_k2, rc_k2) %>% rowMeans(na.rm=T)
  
  avg_m_k2_global <- tmp_df_global %>% 
    dplyr::select(lm_k2, rm_k2) %>% rowMeans(na.rm=T)
  
  avg_c_f_global <- tmp_df_global %>% 
    dplyr::select(lc_f, rc_f) %>% rowMeans(na.rm=T)
  
  avg_m_f_global <- tmp_df_global %>% 
    dplyr::select(lm_f, rm_f) %>% rowMeans(na.rm=T)
  
  avg_c_k2_f_global <- avg_c_k2_global / avg_c_f_global
  avg_m_k2_f_global <- avg_m_k2_global / avg_m_f_global
  
  # Combine results
  results <- bind_cols(
    avg_c_k2_vw, avg_m_k2_vw, avg_c_f_vw, avg_m_f_vw, avg_c_k2_f_vw, avg_m_k2_f_vw,
    avg_c_k2_global, avg_m_k2_global, avg_c_f_global, avg_m_f_global, 
    avg_c_k2_f_global, avg_m_k2_f_global
  ) %>% as.data.frame()
  
  names(results) <- c(
    'cortical_k2_voxel', 'medulla_k2_voxel', 'cortical_f_voxel', 'medulla_f_voxel', 
    'cortical_k2_f_voxel', 'medulla_k2_f_voxel',
    'cortical_k2_global', 'medulla_k2_global', 'cortical_f_global', 'medulla_f_global',
    'cortical_k2_f_global', 'medulla_k2_f_global'
  )
  
  return(results)
}

###All comparison groups 

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

tmp_results <- PET_avg(dat)

dat_results <- dat %>% bind_cols(tmp_results)

dat_results <- dat_results %>% filter(!is.na(cortical_k2_voxel))

dat_results$group2 <- dat_results$group

dat_results$epic_sglti2_1[which(dat_results$group == 'Lean Control')] <- 'No'

RH2 <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/RenalHEIRitage-SGLT2Use.csv')
names(RH2) <- c('Subject', 'event', 'rep_instr', 'rep_inst', 'mrn', 'SGLT2', 'SGLT2_ever')
RH2 <- RH2 %>% filter(!is.na(mrn))

need_med_info <- dat_results %>% filter(group == 'Type 2 Diabetes')
RH2_small <- RH2 %>% filter(mrn %in% need_med_info$mrn)

for(i in c(1:nrow(RH2_small))){
  if(RH2_small$SGLT2[i] == 'No'){
    dat_results$group2[which(dat_results$mrn == RH2_small$mrn[i])] <- 'T2D-No SGLTi2'
    dat_results$epic_sglti2_1[which(dat_results$mrn == RH2_small$mrn[i])] <- 'No'
  }else if(RH2_small$SGLT2[i] == 'Yes'){
    dat_results$group2[which(dat_results$mrn == RH2_small$mrn[i])] <- 'T2D-SGLTi2' 
    dat_results$epic_sglti2_1[which(dat_results$mrn == RH2_small$mrn[i])] <- 'Yes'
  }else{
    next
  }
}

dat_results$epic_sglti2_1[which(dat_results$group2 == 'T2D-SGLTi2')] <- 'Yes'

# Separate test variables for voxel and global
tests_voxel <- c('cortical_k2_voxel', 'cortical_f_voxel', 'cortical_k2_f_voxel',
                 'medulla_k2_voxel', 'medulla_f_voxel', 'medulla_k2_f_voxel')

test_labels_voxel <- c('Cortical K2 (Voxel)', 'Cortical F (Voxel)', 'Cortical K2/F (Voxel)',
                       'Medulla K2 (Voxel)', 'Medulla F (Voxel)', 'Medulla K2/F (Voxel)')

tests_global <- c('cortical_k2_global', 'cortical_f_global', 'cortical_k2_f_global',
                  'medulla_k2_global', 'medulla_f_global', 'medulla_k2_f_global')

test_labels_global <- c('Cortical K2 (Global)', 'Cortical F (Global)', 'Cortical K2/F (Global)',
                        'Medulla K2 (Global)', 'Medulla F (Global)', 'Medulla K2/F (Global)')


boxplot_function <- function(data, variable, label, method){
  
  var_index <- which(names(data) == variable)
  data <- data %>% dplyr::select(group, var_index)
  names(data)[2] <- 'Variable'
  
  desired_order <- c("Lean Control", "Obese Control", "Type 1 Diabetes", "Type 2 Diabetes", "PKD")
  available_levels <- desired_order[desired_order %in% unique(data$group)]
  data$group <- factor(data$group, levels = available_levels)
  
  if(method == 'ANOVA'){
    model <- aov(Variable ~ group, data = data)
    model_results <- TukeyHSD(model, conf.level = 0.95)$group %>% 
      as.data.frame()
    
    model_results <- model_results %>% 
      mutate(pvalue = ifelse(`p adj` < 0.001, '< 0.001', 
                             paste0('p = ', round(`p adj`, 3))))
    
    # Helper function to safely extract p-values
    safe_extract_pval <- function(comparison_name) {
      idx <- which(rownames(model_results) == comparison_name)
      if(length(idx) == 0) {
        return('N/A')
      }
      pval <- model_results$pvalue[idx]
      if(length(pval) == 0 || is.na(pval)) {
        return('N/A')
      }
      return(as.character(pval))
    }
    
    pval_1 <- safe_extract_pval('Type 2 Diabetes-Obese Control')
    pval_2 <- safe_extract_pval('Type 2 Diabetes-Lean Control')
    pval_3 <- safe_extract_pval('Type 2 Diabetes-Type 1 Diabetes')
    pval_4 <- safe_extract_pval('PKD-Type 2 Diabetes')
    
  }else if(method == 't-test'){
    
    # Helper function to perform t-test with error handling
    safe_ttest <- function(data, groups) {
      tmp_df <- data %>% filter(group %in% groups, !is.na(Variable))
      
      # Check if both groups exist and have data
      group_counts <- table(tmp_df$group)
      if(length(group_counts) != 2 || any(group_counts < 2)) {
        return('N/A')
      }
      
      tryCatch({
        model1 <- t.test(Variable ~ group, data = tmp_df)
        pval <- ifelse(model1$p.value < 0.001, '< 0.001', 
                       paste0('p = ', round(model1$p.value, 3)))
        return(pval)
      }, error = function(e) {
        return('N/A')
      })
    }
    
    pval_1 <- safe_ttest(data, c('Type 2 Diabetes', 'Obese Control'))
    pval_2 <- safe_ttest(data, c('Type 2 Diabetes', 'Lean Control'))
    pval_3 <- safe_ttest(data, c('Type 2 Diabetes', 'Type 1 Diabetes'))
    pval_4 <- safe_ttest(data, c('Type 2 Diabetes', 'PKD'))
  }
  
  y_max <- max(data$Variable, na.rm = TRUE)
  y_range <- diff(range(data$Variable, na.rm = TRUE))
  
  plot <- ggplot(data, aes(x = group, y = Variable, fill = group))  +
    geom_boxplot()+
    scale_fill_manual(values = c("Lean Control" = "#87CEEB", 
                                 "Obese Control" = "#ADD8E6", 
                                 "Type 1 Diabetes" = "#F0E68C", 
                                 "Type 2 Diabetes" = "#CD5C5C", 
                                 "PKD" = "#DDA0DD")) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) + 
    labs(x= 'Study Group', y = label, fill = 'Study Group')+
    theme_minimal()+
    theme(axis.text.x = element_blank(),
          text = element_text(size = 20),
          legend.position = "right",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  group_positions <- 1:length(levels(data$group))
  names(group_positions) <- levels(data$group)
  
  get_position <- function(group_name) {
    pos <- group_positions[group_name]
    if(is.na(pos)) {
      return(NULL)
    }
    return(as.numeric(pos))
  }
  
  t2d_pos <- get_position("Type 2 Diabetes")
  lean_pos <- get_position("Lean Control")
  obese_pos <- get_position("Obese Control")
  t1d_pos <- get_position("Type 1 Diabetes")
  pkd_pos <- get_position("PKD")
  
  # Helper function to check if we should add comparison
  should_add_comparison <- function(pos1, pos2, pval) {
    !is.null(pos1) && !is.null(pos2) && !is.na(pval) && pval != 'N/A'
  }
  
  # Only add comparison lines if positions exist and p-values are available
  if(should_add_comparison(lean_pos, t2d_pos, pval_2)){
    plot <- plot + 
      annotate("segment", x = lean_pos, xend = t2d_pos, 
               y = y_max + 0.25 * y_range, yend = y_max + 0.25 * y_range, 
               color = "black", size = 0.5) +
      annotate("segment", x = lean_pos, xend = lean_pos, 
               y = y_max + 0.23 * y_range, yend = y_max + 0.25 * y_range, 
               color = "black", size = 0.5) +
      annotate("segment", x = t2d_pos, xend = t2d_pos, 
               y = y_max + 0.23 * y_range, yend = y_max + 0.25 * y_range, 
               color = "black", size = 0.5) +
      annotate("text", x = (lean_pos + t2d_pos)/2, y = y_max + 0.27 * y_range, 
               label = pval_2, size = 4)
  }
  
  if(should_add_comparison(obese_pos, t2d_pos, pval_1)){
    plot <- plot +
      annotate("segment", x = obese_pos, xend = t2d_pos, 
               y = y_max + 0.18 * y_range, yend = y_max + 0.18 * y_range, 
               color = "black", size = 0.5) +
      annotate("segment", x = obese_pos, xend = obese_pos, 
               y = y_max + 0.16 * y_range, yend = y_max + 0.18 * y_range, 
               color = "black", size = 0.5) +
      annotate("segment", x = t2d_pos, xend = t2d_pos, 
               y = y_max + 0.16 * y_range, yend = y_max + 0.18 * y_range, 
               color = "black", size = 0.5) +
      annotate("text", x = (obese_pos + t2d_pos)/2, y = y_max + 0.20 * y_range, 
               label = pval_1, size = 4)
  }
  
  if(should_add_comparison(t1d_pos, t2d_pos, pval_3)){
    plot <- plot +
      annotate("segment", x = t1d_pos, xend = t2d_pos, 
               y = y_max + 0.11 * y_range, yend = y_max + 0.11 * y_range, 
               color = "black", size = 0.5) +
      annotate("segment", x = t1d_pos, xend = t1d_pos, 
               y = y_max + 0.09 * y_range, yend = y_max + 0.11 * y_range, 
               color = "black", size = 0.5) +
      annotate("segment", x = t2d_pos, xend = t2d_pos, 
               y = y_max + 0.09 * y_range, yend = y_max + 0.11 * y_range, 
               color = "black", size = 0.5) +
      annotate("text", x = (t1d_pos + t2d_pos)/2, y = y_max + 0.13 * y_range, 
               label = pval_3, size = 4)
  }
  
  if(should_add_comparison(t2d_pos, pkd_pos, pval_4)){
    plot <- plot +
      annotate("segment", x = t2d_pos, xend = pkd_pos, 
               y = y_max + 0.04 * y_range, yend = y_max + 0.04 * y_range, 
               color = "black", size = 0.5) +
      annotate("segment", x = t2d_pos, xend = t2d_pos, 
               y = y_max + 0.02 * y_range, yend = y_max + 0.04 * y_range, 
               color = "black", size = 0.5) +
      annotate("segment", x = pkd_pos, xend = pkd_pos, 
               y = y_max + 0.02 * y_range, yend = y_max + 0.04 * y_range, 
               color = "black", size = 0.5) +
      annotate("text", x = (t2d_pos + pkd_pos)/2, y = y_max + 0.06 * y_range, 
               label = pval_4, size = 4)
  }
  
  plot <- plot + expand_limits(y = y_max + 0.3 * y_range)
  
  return(plot)
}

# ============== VOXEL PLOTS - ANOVA ==============
for(i in c(1:length(tests_voxel))){
  if(i == 1){
    results_list_voxel <- list()
  }
  results_list_voxel[[i]] <- boxplot_function(dat_results, tests_voxel[i], 
                                              test_labels_voxel[i], method='ANOVA')
}

png('C:/Users/netio/Documents/UofW/Projects/Imaging_Shivani/AllComparisonGroups_KidneyImaging_VOXEL.png', 
    width = 1200, height = 1600)  
gridExtra::grid.arrange(grobs = results_list_voxel, ncol = 2)
dev.off()

# ============== VOXEL PLOTS - T-TEST ==============
for(i in c(1:length(tests_voxel))){
  if(i == 1){
    results_list_voxel <- list()
  }
  results_list_voxel[[i]] <- boxplot_function(dat_results, tests_voxel[i], 
                                              test_labels_voxel[i], method='t-test')
}

png('C:/Users/netio/Documents/UofW/Projects/Imaging_Shivani/AllComparisonGroups_KidneyImaging_VOXEL_ttest.png', 
    width = 1200, height = 1600)  
gridExtra::grid.arrange(grobs = results_list_voxel, ncol = 2)
dev.off()

# ============== GLOBAL PLOTS - ANOVA ==============
for(i in c(1:length(tests_global))){
  if(i == 1){
    results_list_global <- list()
  }
  results_list_global[[i]] <- boxplot_function(dat_results, tests_global[i], 
                                               test_labels_global[i], method='ANOVA')
}

png('C:/Users/netio/Documents/UofW/Projects/Imaging_Shivani/AllComparisonGroups_KidneyImaging_GLOBAL.png', 
    width = 1200, height = 1600)  
gridExtra::grid.arrange(grobs = results_list_global, ncol = 2)
dev.off()

# ============== GLOBAL PLOTS - T-TEST ==============
for(i in c(1:length(tests_global))){
  if(i == 1){
    results_list_global <- list()
  }
  results_list_global[[i]] <- boxplot_function(dat_results, tests_global[i], 
                                               test_labels_global[i], method='t-test')
}

png('C:/Users/netio/Documents/UofW/Projects/Imaging_Shivani/AllComparisonGroups_KidneyImaging_GLOBAL_ttest.png', 
    width = 1200, height = 1600)  
gridExtra::grid.arrange(grobs = results_list_global, ncol = 2)
dev.off()








