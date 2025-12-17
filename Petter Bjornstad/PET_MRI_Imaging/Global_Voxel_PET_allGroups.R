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
  names(results) <- c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 
                      'avg_c_k2_f', 'avg_m_k2_f')
  
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


dat_results <- dat


dat_results <- dat_results %>% bind_cols(tmp_results)


dat_results <- dat_results %>% filter(!is.na(avg_c_k2))

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





tests <- c('avg_c_k2', 'avg_c_f', 
           'avg_m_k2', 'avg_m_f', 
           'avg_c_k2_f', 'avg_m_k2_f')



boxplot_function <- function(data, variable, label, method){
  
  var_index <- which(names(data) == variable)
  data <- data %>% dplyr::select(group, var_index)
  names(data)[2] <- 'Variable'
  
  
  desired_order <- c("Lean Control", "Obese Control", "Type 1 Diabetes", "Type 2 Diabetes", "PKD")
  
  # Only include levels that actually exist in your data
  available_levels <- desired_order[desired_order %in% unique(data$group)]
  data$group <- factor(data$group, levels = available_levels)
  
  if(method == 'ANOVA'){
    model <- aov(Variable ~ group, data = data)
    model_results <- TukeyHSD(model, conf.level = 0.95)$group %>% 
      as.data.frame()
    
    model_results <- model_results %>% 
      mutate(pvalue = ifelse(`p adj` < 0.001, '< 0.001', 
                             paste0('p = ', round(`p adj`, 3))))
    
    pval_1 <- model_results$pvalue[which(rownames(model_results) == 'Type 2 Diabetes-Obese Control')] %>%
      as.character()
    pval_2 <- model_results$pvalue[which(rownames(model_results) == 'Type 2 Diabetes-Lean Control')] %>% 
      as.character()
    pval_3 <- model_results$pvalue[which(rownames(model_results) == 'Type 2 Diabetes-Type 1 Diabetes')] %>% 
      as.character()
    pval_4 <- model_results$pvalue[which(rownames(model_results) == 'PKD-Type 2 Diabetes')] %>% 
      as.character()
  }else if(method == 't-test'){
    
    tmp_df <- data %>% filter(group %in% c('Type 2 Diabetes', 'Obese Control'))
    model1 <- t.test(Variable ~ group, data = tmp_df)
    pval_1 <- ifelse(model1$p.value < 0.001, '< 0.001', 
                     paste0('p = ', round(model1$p.value, 3)))
    
    tmp_df <- data %>% filter(group %in% c('Type 2 Diabetes', 'Lean Control'))
    model1 <- t.test(Variable ~ group, data = tmp_df)
    pval_2 <- ifelse(model1$p.value < 0.001, '< 0.001', 
                     paste0('p = ', round(model1$p.value, 3)))
    
    tmp_df <- data %>% filter(group %in% c('Type 2 Diabetes', 'Type 1 Diabetes'))
    model1 <- t.test(Variable ~ group, data = tmp_df)
    pval_3 <- ifelse(model1$p.value < 0.001, '< 0.001', 
                     paste0('p = ', round(model1$p.value, 3)))
    
    tmp_df <- data %>% filter(group %in% c('Type 2 Diabetes', 'PKD'))
    model1 <- t.test(Variable ~ group, data = tmp_df)
    pval_4 <- ifelse(model1$p.value < 0.001, '< 0.001', 
                     paste0('p = ', round(model1$p.value, 3)))
    
    
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
    theme(axis.text.x = element_blank(),  # Keep x-axis labels hidden as in your original
          text = element_text(size = 20),
          legend.position = "right",      # Show legend since x-axis is hidden
          panel.grid.major.x = element_blank(),  # Remove vertical grid lines
          panel.grid.minor.x = element_blank())
  
  group_positions <- 1:length(levels(data$group))
  names(group_positions) <- levels(data$group)
  
  # Debug: Print positions
  print("Group positions:")
  print(group_positions)
  
  # Find positions of specific groups (with error checking)
  get_position <- function(group_name) {
    pos <- group_positions[group_name]
    if(is.na(pos)) {
      print(paste("Warning: Group", group_name, "not found in data"))
      return(NULL)
    }
    return(as.numeric(pos))
  }
  
  t2d_pos <- get_position("Type 2 Diabetes")
  lean_pos <- get_position("Lean Control")
  obese_pos <- get_position("Obese Control")
  t1d_pos <- get_position("Type 1 Diabetes")
  pkd_pos <- get_position("PKD")
  
  
  
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
             label = pval_2, size = 4) +
    
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
             label = pval_1, size = 4) +
    
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
             label = pval_3, size = 4) +
    
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
             label = pval_4, size = 4) + 
    
    expand_limits(y = y_max + 0.3 * y_range)
  
  print(plot)
}






for(i in c(1:length(tests))){
  if(i == 1){
    results_list <- list()
  }
  results_list[[i]]<- boxplot_function(dat_results, tests[i], tests[i], method='ANOVA')
}

png('C:/Users/netio/Documents/UofW/Rockies/AllComparisonGroups_KidneyImaging.png', 
    width =1200, height = 1600)  
gridExtra::grid.arrange(results_list[[1]], results_list[[2]], 
                        results_list[[3]], results_list[[4]], 
                        results_list[[5]], results_list[[6]],
                        ncol = 2)

dev.off()



for(i in c(1:length(tests))){
  if(i == 1){
    results_list <- list()
  }
  results_list[[i]]<- boxplot_function(dat_results, tests[i], tests[i], method='t-test')
}

png('C:/Users/netio/Documents/UofW/Rockies/AllComparisonGroups_KidneyImaging_ttest.png', 
    width =1200, height = 1600)  
gridExtra::grid.arrange(results_list[[1]], results_list[[2]], 
                        results_list[[3]], results_list[[4]], 
                        results_list[[5]], results_list[[6]],
                        ncol = 2)

dev.off()









