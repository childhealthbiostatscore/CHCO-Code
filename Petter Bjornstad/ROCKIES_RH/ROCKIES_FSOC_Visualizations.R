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





load('C:/Users/netio/Documents/UofW/Rockies/Line438_Boxplots_NoMed.RData')


harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))



dat2 <- dat %>% filter(visit == 'baseline') %>% 
#  filter(study %in% c('RENAL-HEIR', 'RENAL-HEIRitage', 'CROCODILE') | record_id == 'IT_19') %>%
#  filter(group != 'Obese Control') %>% 
  dplyr::select(mrn, record_id, study, visit, group, starts_with('fsoc'), bmi, hba1c,
                epic_sglti2_1, age, sex, epic_mfm_1, epic_insulin_1, epic_glp1ra_1)

dat2 <- dat2[-which(dat2$group == 'Type 1 Diabetes'),]


tests <- c('fsoc_l_cortex', 'fsoc_r_cortex', 
           'fsoc_l_kidney', 'fsoc_r_kidney', 
           'fsoc_l_medulla', 'fsoc_r_medulla', 
           'fsoc_medulla', 'fsoc_cortex', 'fsoc_kidney')


graphs <- list()


dat2 <- dat2 %>% filter(!is.na(fsoc_l_cortex) | !is.na(fsoc_r_cortex))

dat2$group2 <- dat2$group

dat2$group2[which(dat2$group == 'Type 2 Diabetes' & dat2$epic_sglti2_1 == 'Yes')] <- 'T2D-SGLTi2'
dat2$group2[which(dat2$group == 'Type 2 Diabetes' & dat2$epic_sglti2_1 == 'No')] <- 'T2D-No SGLTi2'
dat2$epic_sglti2_1[which(dat2$group == 'Lean Control')] <- 'No'
dat2$epic_sglti2_1[which(dat2$group == 'Obese Control')] <- 'No'


medications <- readxl::read_xlsx("C:/Users/netio/Documents/Harmonized_data/Biopsies_w_mrn_Oct3.xlsx")
medications <- medications %>% dplyr::select(mrn, ends_with('_1'), -starts_with('ever_'))
names(medications) <- str_replace(names(medications), pattern = '_1', replacement = '')



need_med_info <- dat2 %>% filter(group == 'Type 2 Diabetes') %>% filter(is.na(epic_sglti2_1))

med_small <- medications %>% filter(mrn %in% need_med_info$mrn)

for(i in c(1:nrow(med_small))){
  if(med_small$sglti2[i] == 'No'){
    dat2$group2[which(dat2$mrn == med_small$mrn[i])] <- 'T2D-No SGLTi2'
    dat2$epic_sglti2_1[which(dat2$mrn == med_small$mrn[i])] <- 'No'
  }else if(med_small$sglti2[i] == 'Yes'){
    dat2$group2[which(dat2$mrn == med_small$mrn[i])] <- 'T2D-SGLTi2'
    dat2$epic_sglti2_1[which(dat2$mrn == med_small$mrn[i])] <- 'No'
  }else{
    next
  }
}

med_small$mrn <- as.numeric(med_small$mrn)
need_med_info <- need_med_info %>% anti_join(med_small, by='mrn')

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

dat2$group2[which(dat2$group == 'Obese Control')] <- 'Obese Control' 
dat2$group2[which(dat2$group == 'Lean Control')] <- 'Lean Control' 



dat2 <- dat2 %>% filter(!is.na(fsoc_l_cortex) | !is.na(fsoc_r_cortex))



tmp_df <- dat2 %>% dplyr::select(starts_with('fsoc'))



find_fsoc_averages <- function(data){
  tmp_data <- data %>% dplyr::select(starts_with('fsoc'))
  fsoc_full_combined <- rowMeans(tmp_data, na.rm=T)
  
  tmp_data <- data %>% dplyr::select(starts_with('fsoc_l_'))
  fsoc_l_combined <- tmp_data %>% rowMeans(na.rm=T)
  
  tmp_data <- data %>% dplyr::select(starts_with('fsoc_r_'))
  fsoc_r_combined <- tmp_data %>% rowMeans(na.rm=T)
  
  fsoc_medulla <- data %>% dplyr::select(fsoc_l_medulla, fsoc_r_medulla) %>%
    rowMeans(na.rm=T)
  fsoc_cortex <- data %>% dplyr::select(fsoc_l_medulla, fsoc_r_cortex) %>%
    rowMeans(na.rm=T)
  fsoc_kidney <- data %>% dplyr::select(fsoc_l_medulla, fsoc_r_kidney) %>%
    rowMeans(na.rm=T)
  
  tmp_df <- cbind(fsoc_l_combined, fsoc_r_combined, fsoc_medulla, fsoc_cortex, fsoc_kidney, fsoc_full_combined)
  return(tmp_df)
  
}

tmp_results <- find_fsoc_averages(dat2)


dat_results <- dat2 %>% bind_cols(tmp_results)









apply(dat_results %>% dplyr::select(starts_with('fsoc')), 2, is.na) %>% 
  colSums() %>% which.min()



#pdf('C:/Users/netio/Documents/UofW/Rockies/SGLT2ComparisonGroups_KidneyImaging.pdf', 
#    width =15, height = 15)

dat_results_combined <- dat_results %>% 
  filter(group2 %in% c('T2D-No SGLTi2', 'T2D-SGLTi2')) %>% 
  mutate(group2 = 'T2D Combined')

df_plot <- bind_rows(dat_results, dat_results_combined) %>% 
  mutate(group2 = factor(group2, levels = c('Lean Control',
                                            'Obese Control',
                                            'T2D-No SGLTi2', 
                                            'T2D-SGLTi2',
                                            'T2D Combined')))


df_plot <- df_plot %>% 
  dplyr::select(mrn, group, group2, fsoc_l_cortex, fsoc_r_cortex, 
                fsoc_l_kidney, fsoc_r_kidney, fsoc_l_medulla, fsoc_r_medulla,
                fsoc_l_combined, fsoc_r_combined, fsoc_full_combined, 
                fsoc_medulla, fsoc_cortex, fsoc_kidney)


df_plot <- df_plot %>% mutate(position = ifelse(group2 == 'T2D-No SGLTi2', 1.7, ifelse(group2 == 'T2D-SGLTi2', 2.0, NA)))



table1::table1(~age + sex + bmi+ hba1c + study + epic_mfm_1+epic_insulin_1+epic_glp1ra_1+epic_sglti2_1 + fsoc_l_kidney + fsoc_r_kidney | group2,
               data=dat_results %>% 
                 filter(!is.na(fsoc_l_cortex) | !is.na(fsoc_r_cortex)))



boxplot_function <- function(data, variable, label, method){
  
  var_index <- which(names(data) == variable)
  data <- data %>% dplyr::select(group2, var_index)
  names(data)[2] <- 'Variable'
  
  data <- data %>% 
    mutate(position = case_when(
      group2 == 'Lean Control' ~ 1,
      group2 == 'Obese Control' ~ 2,
      group2 == 'T2D-No SGLTi2' ~ 2.8,
      group2 == 'T2D-SGLTi2' ~ 3.2,
      group2 == 'T2D Combined' ~ 3, 
      TRUE ~ NA_real_))
  
  if(method == 'ANOVA'){
  model <- aov(Variable ~ group2, data = data)
  model_results <- TukeyHSD(model, conf.level = 0.95)$group2 %>% 
    as.data.frame()
  
  model_results <- model_results %>% 
    mutate(pvalue = ifelse(`p adj` < 0.001, '< 0.001', 
                           paste0('p = ', round(`p adj`, 3))))
  
  pval_T2D_noslgt2_lean <- model_results$pvalue[which(rownames(model_results) == 'T2D-No SGLTi2-Lean Control')] %>%
    as.character()
  pval_T2D_total_lean <- model_results$pvalue[which(rownames(model_results) == 'T2D Combined-Lean Control')] %>% 
    as.character()
  pval_T2D_noslgt2_obese <- model_results$pvalue[which(rownames(model_results) == 'T2D-No SGLTi2-Obese Control')] %>%
    as.character()
  pval_T2D_total_obese <- model_results$pvalue[which(rownames(model_results) == 'T2D Combined-Obese Control')] %>% 
    as.character()
  pval_T2D_comparison <- model_results$pvalue[which(rownames(model_results) == 'T2D-SGLTi2-T2D-No SGLTi2')] %>% 
    as.character()
  
  }else if(method == 't-test'){
    tmp_df <- data %>% filter(group2 %in% c('T2D-No SGLTi2', 'Lean Control'))
    model1 <- t.test(Variable ~ group2, data = tmp_df)
    pval_T2D_noslgt2_lean <- ifelse(model1$p.value < 0.001, '< 0.001', 
                                       paste0('p = ', round(model1$p.value, 3)))
    
    tmp_df <- data %>% filter(group2 %in% c('T2D Combined', 'Lean Control'))
    model1 <- t.test(Variable ~ group2, data = tmp_df)
    pval_T2D_total_lean <- ifelse(model1$p.value < 0.001, '< 0.001', 
                                       paste0('p = ', round(model1$p.value, 3)))
    
    tmp_df <- data %>% filter(group2 %in% c('T2D-No SGLTi2', 'Obese Control'))
    model1 <- t.test(Variable ~ group2, data = tmp_df)
    pval_T2D_noslgt2_obese <- ifelse(model1$p.value < 0.001, '< 0.001', 
                                    paste0('p = ', round(model1$p.value, 3)))
    
    tmp_df <- data %>% filter(group2 %in% c('T2D Combined', 'Obese Control'))
    model1 <- t.test(Variable ~ group2, data = tmp_df)
    pval_T2D_total_obese <- ifelse(model1$p.value < 0.001, '< 0.001', 
                                  paste0('p = ', round(model1$p.value, 3)))
    
    tmp_df <- data %>% filter(group2 %in% c('T2D-No SGLTi2', 'T2D-SGLTi2'))
    model1 <- t.test(Variable ~ group2, data = tmp_df)
    pval_T2D_comparison <- ifelse(model1$p.value < 0.001, '< 0.001', 
                                     paste0('p = ', round(model1$p.value, 3)))
    
    
    
  }
  
  y_max <- max(data$Variable, na.rm = TRUE)
  y_range <- diff(range(data$Variable, na.rm = TRUE))

    
  plot <- ggplot(data %>% dplyr::filter(group2 %in% c('Lean Control', 'Obese Control', 'T2D Combined')), 
                 aes(x = position, y = Variable, fill = group2))  +
    geom_boxplot(width = 0.6, size = 1)+
    scale_fill_manual(values = c("#c2dfe3", "#9bc1bc", "#fff9ec", "#fcb1a6", "#fb6376")) +
    geom_boxplot(data = data %>%
                   dplyr::filter(group2 %in% c('T2D-No SGLTi2', 'T2D-SGLTi2')), 
                 aes(x = position, y=Variable, fill = group2), width = 0.1, size = 1)+
    labs(x= 'Study Group', y = label, fill = 'Study Group')+
    scale_x_continuous(breaks = c(1, 2, 2.8, 3, 3.2), 
                       labels = c('Lean Control', 'Obese Control', 'T2D-No SGLT2', 'T2D Combined', 'T2D-SGLT2'),
                       limits = c(0.5, 4.1), 
                       expand = expansion(mult = c(0.05, 0.05)))+
    theme_minimal()+
    theme(axis.text.x = element_blank(),
          text = element_text(size = 20))
  
  plot <- plot + 
    # T2D Combined vs Lean Control (top level)
    geom_segment(aes(x = 1, xend = 3, y = y_max + 0.25 * y_range, yend = y_max + 0.25 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 1, xend = 1, y = y_max + 0.23 * y_range, yend = y_max + 0.25 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 3, xend = 3, y = y_max + 0.23 * y_range, yend = y_max + 0.25 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    annotate("text", x = 2, y = y_max + 0.29 * y_range, label = pval_T2D_total_lean, size = 4.5) +
    
    # T2D Combined vs Obese Control (second level)
    geom_segment(aes(x = 2, xend = 3, y = y_max + 0.19 * y_range, yend = y_max + 0.19 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 2, xend = 2, y = y_max + 0.17 * y_range, yend = y_max + 0.19 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 3, xend = 3, y = y_max + 0.17 * y_range, yend = y_max + 0.19 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    annotate("text", x = 2.5, y = y_max + 0.23 * y_range, label = pval_T2D_total_obese, size = 4.5) +
    
    # T2D-No SGLTi2 vs Lean Control
    geom_segment(aes(x = 1, xend = 2.7, y = y_max + 0.13 * y_range, yend = y_max + 0.13 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 1, xend = 1, y = y_max + 0.11 * y_range, yend = y_max + 0.13 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 2.7, xend = 2.7, y = y_max + 0.11 * y_range, yend = y_max + 0.13 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    annotate("text", x = 1.85, y = y_max + 0.17 * y_range, label = pval_T2D_noslgt2_lean, size = 4.5) +
    
    # T2D-No SGLTi2 vs Obese Control
    geom_segment(aes(x = 2, xend = 2.7, y = y_max + 0.07 * y_range, yend = y_max + 0.07 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 2, xend = 2, y = y_max + 0.05 * y_range, yend = y_max + 0.07 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 2.7, xend = 2.7, y = y_max + 0.05 * y_range, yend = y_max + 0.07 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    annotate("text", x = 2.35, y = y_max + 0.11 * y_range, label = pval_T2D_noslgt2_obese, size = 4.5) +
    
    # T2D subgroup comparison (bottom level)
    geom_segment(aes(x = 2.8, xend = 3.2, y = y_max + 0.01 * y_range, yend = y_max + 0.01 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 2.8, xend = 2.8, y = y_max - 0.01 * y_range, yend = y_max + 0.01 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 3.2, xend = 3.2, y = y_max - 0.01 * y_range, yend = y_max + 0.01 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    annotate("text", x = 3.25, y = y_max + 0.05 * y_range, label = pval_T2D_comparison, size = 4.5) +
    expand_limits(y = y_max + 0.35 * y_range)
  
  print(plot)
}




for(i in c(1:length(tests))){
  if(i == 1){
    results_list <- list()
  }
  results_list[[i]]<- boxplot_function(df_plot, tests[i], tests[i], method='ANOVA')
}


library(grid)


pdf('C:/Users/netio/Documents/UofW/Rockies/SGLT2ComparisonGroups_FSOC.pdf', 
    width =20, height = 20)  
gridExtra::grid.arrange(results_list[[1]], results_list[[2]], 
                        results_list[[3]], results_list[[4]], 
                        results_list[[5]], results_list[[6]], 
                        results_list[[7]], results_list[[8]],
                        results_list[[9]],
                        ncol = 2)

dev.off()

png('C:/Users/netio/Documents/UofW/Rockies/SGLT2ComparisonGroups_FSOC.png', 
    width =1200, height = 1600)  
gridExtra::grid.arrange(results_list[[1]], results_list[[2]], 
                        results_list[[3]], results_list[[4]], 
                        results_list[[5]], results_list[[6]], 
                        results_list[[7]], results_list[[8]],
                        results_list[[9]],
                        ncol = 2)

dev.off()




#t-tests analysis 



for(i in c(1:length(tests))){
  if(i == 1){
    results_list <- list()
  }
  results_list[[i]]<- boxplot_function(df_plot, tests[i], tests[i], method='t-test')
}


library(grid)


pdf('C:/Users/netio/Documents/UofW/Rockies/SGLT2ComparisonGroups_FSOC_ttest.pdf', 
    width =20, height = 20)  
gridExtra::grid.arrange(results_list[[1]], results_list[[2]], 
                        results_list[[3]], results_list[[4]], 
                        results_list[[5]], results_list[[6]], 
                        results_list[[7]], results_list[[8]],
                        results_list[[9]],
                        ncol = 2)

dev.off()

png('C:/Users/netio/Documents/UofW/Rockies/SGLT2ComparisonGroups_FSOC_ttest.png', 
    width =1200, height = 1600)  
gridExtra::grid.arrange(results_list[[1]], results_list[[2]], 
                        results_list[[3]], results_list[[4]], 
                        results_list[[5]], results_list[[6]], 
                        results_list[[7]], results_list[[8]],
                        results_list[[9]],
                        ncol = 2)

dev.off()













