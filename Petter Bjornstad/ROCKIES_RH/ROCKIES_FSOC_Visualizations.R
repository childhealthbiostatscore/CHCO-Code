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

dat2 <- dat %>% filter(visit == 'baseline') %>% 
  filter(study %in% c('RENAL-HEIR', 'RENAL-HEIRitage', 'CROCODILE') | record_id == 'IT_19') %>%
#  filter(group != 'Obese Control') %>% 
  dplyr::select(mrn, record_id, study, visit, group, group2, starts_with('fsoc'), bmi, 
                epic_sglti2_1, age, sex, epic_mfm_1, epic_insulin_1, epic_glp1ra_1)

dat2 <- dat2[-which(dat2$group == 'Type 1 Diabetes'),]


tests <- c('fsoc_l_cortex', 'fsoc_r_cortex', 
           'fsoc_l_kidney', 'fsoc_r_kidney', 
           'fsoc_l_medulla', 'fsoc_r_medulla', 
           'fsoc_l_combined', 'fsoc_r_combined',
           'fsoc_medulla', 'fsoc_cortex', 'fsoc_kidney',
           'fsoc_full_combined')


graphs <- list()




medications <- readxl::read_xlsx("C:/Users/netio/Documents/Harmonized_data/Biopsies_w_mrn_Oct3.xlsx")
medications <- medications %>% dplyr::select(mrn, ends_with('_1'), -starts_with('ever_'))
names(medications) <- str_replace(names(medications), pattern = '_1', replacement = '')



need_med_info <- dat2 %>% filter(group == 'Type 2 Diabetes') %>% filter(is.na(group2))

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

dat2$epic_sglti2_1[which(dat2$group == 'Lean Control')] <- 'No'

dat2$group2[which(dat2$group == 'Obese Control')] <- 'Obese Control' 


table1::table1(~age + sex + bmi+study + epic_mfm_1+epic_insulin_1+epic_glp1ra_1+epic_sglti2_1 | group2,data=dat2 %>% 
                 filter(!is.na(group2)))

dat2$group2[which(dat2$group %in% c('Obese Control', 'Lean Control'))] <- 'Control'


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



max_num_trait <- apply(dat_results %>% dplyr::select(starts_with('fsoc')), 2, is.na) %>% colSums() %>% which.min() %>% names()

min_num_missing <- sum(is.na(dplyr::select(dat_results, starts_with('fsoc'))[,apply(dat_results %>% dplyr::select(starts_with('fsoc')), 2, is.na) %>% colSums() %>% which.min()]))


max_num_lc <- dat_results %>% filter(group == 'Lean Control') %>% filter(!is.na(fsoc_r_cortex)) %>% 
  nrow()
max_num_nosglt2 <- dat_results %>% filter(group == 'T2D-No SGLTi2') %>% filter(!is.na(fsoc_r_cortex)) %>% 
  nrow()
max_num_sglt2 <- dat_results %>% filter(group == 'T2D-SGLTi2') %>% filter(!is.na(fsoc_r_cortex)) %>% 
  nrow()



number_part_df <- dat_results %>% filter(!is.na(fsoc_full_combined))

min_num_lc <- number_part_df %>% filter(group2 == 'Lean Control')
min_num_nosglt2 <- number_part_df %>% filter(group2 == 'T2D-No SGLTi2')
min_num_sglt2 <- number_part_df %>% filter(group2 == 'T2D-SGLTi2')




#pdf('C:/Users/netio/Documents/UofW/Rockies/SGLT2ComparisonGroups_KidneyImaging.pdf', 
#    width =15, height = 15)

dat_results_combined <- dat_results %>% 
  filter(group2 %in% c('T2D-No SGLTi2', 'T2D-SGLTi2')) %>% 
  mutate(group2 = 'T2D Combined')

df_plot <- bind_rows(dat_results, dat_results_combined) %>% 
  mutate(group2 = factor(group2, levels = c('Control', 
                                            'T2D-No SGLTi2', 
                                            'T2D-SGLTi2',
                                            'T2D Combined')))


df_plot <- df_plot %>% 
  dplyr::select(mrn, group, group2, fsoc_l_cortex, fsoc_r_cortex, 
                fsoc_l_kidney, fsoc_r_kidney, fsoc_l_medulla, fsoc_r_medulla,
                fsoc_l_combined, fsoc_r_combined, fsoc_full_combined, 
                fsoc_medulla, fsoc_cortex, fsoc_kidney)


df_plot <- df_plot %>% mutate(position = ifelse(group2 == 'T2D-No SGLTi2', 1.7, ifelse(group2 == 'T2D-SGLTi2', 2.0, NA)))



table1::table1(~age + sex + bmi+study + epic_mfm_1+epic_insulin_1+epic_glp1ra_1+epic_sglti2_1 + fsoc_l_cortex | group2,data=dat_results %>% 
                 filter(!is.na(fsoc_l_cortex)))



boxplot_function <- function(data, variable, label, method){
  
  var_index <- which(names(data) == variable)
  data <- data %>% dplyr::select(group2, var_index)
  names(data)[2] <- 'Variable'
  
  data <- data %>% mutate(position = ifelse(group2 == 'T2D-No SGLTi2', 1.7, ifelse(group2 == 'T2D-SGLTi2', 2.0, NA)))
  
  if(method == 'ANOVA'){
  model <- aov(Variable ~ group2, data = data)
  model_results <- TukeyHSD(model, conf.level = 0.95)$group2 %>% 
    as.data.frame()
  
  model_results <- model_results %>% 
    mutate(pvalue = ifelse(`p adj` < 0.001, '< 0.001', 
                           paste0('p = ', round(`p adj`, 3))))
  
  pval_T2D_noslgt2_control <- model_results$pvalue[which(rownames(model_results) == 'T2D-No SGLTi2-Control')] %>%
    as.character()
  pval_T2D_total_control <- model_results$pvalue[which(rownames(model_results) == 'T2D Combined-Control')] %>% 
    as.character()
  pval_T2D_comparison <- model_results$pvalue[which(rownames(model_results) == 'T2D-SGLTi2-T2D-No SGLTi2')] %>% 
    as.character()
  }else if(method == 't-test'){
    tmp_df <- data %>% filter(group2 %in% c('T2D-No SGLTi2', 'Control'))
    model1 <- t.test(Variable ~ group2, data = tmp_df)
    pval_T2D_noslgt2_control <- ifelse(model1$p.value < 0.001, '< 0.001', 
                                       paste0('p = ', round(model1$p.value, 3)))
    
    tmp_df <- data %>% filter(group2 %in% c('T2D Combined', 'Control'))
    model1 <- t.test(Variable ~ group2, data = tmp_df)
    pval_T2D_total_control <- ifelse(model1$p.value < 0.001, '< 0.001', 
                                       paste0('p = ', round(model1$p.value, 3)))
    
    tmp_df <- data %>% filter(group2 %in% c('T2D-No SGLTi2', 'T2D-SGLTi2'))
    model1 <- t.test(Variable ~ group2, data = tmp_df)
    pval_T2D_comparison <- ifelse(model1$p.value < 0.001, '< 0.001', 
                                     paste0('p = ', round(model1$p.value, 3)))
    
    
    
  }
  
  y_max <- max(data$Variable, na.rm = TRUE)
  y_range <- diff(range(data$Variable, na.rm = TRUE))

    
  plot <- ggplot(data %>% dplyr::filter(group2 %in% c('Control', 'T2D Combined')), aes(x = group2, y = Variable, fill = group2))  +
    geom_boxplot(width = 1.3, size = 1)+
    scale_fill_manual(values = c("#c2dfe3", "#fff9ec", "#fcb1a6", "#fb6376")) +
    geom_boxplot(data = data %>%
                   dplyr::filter(group2 %in% c('T2D-No SGLTi2', 'T2D-SGLTi2')), 
                 aes(x = position, y=Variable, fill = group2), width = 0.1, size = 1)+
    labs(x= 'Study Group', y = label, fill = 'Study Group')+
    theme_minimal()+
    theme(axis.text.x = element_blank(),
          text = element_text(size = 20))
  
  plot <- plot + 
    geom_segment(aes(x = 1, xend = 2, y = y_max + 0.15 * y_range, yend = y_max + 0.15 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 1, xend = 1, y = y_max + 0.13 * y_range, yend = y_max + 0.15 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 2, xend = 2, y = y_max + 0.13 * y_range, yend = y_max + 0.15 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    annotate("text", x = 1.5, y = y_max + 0.19 * y_range, label = pval_T2D_total_control, size = 4.5) +
    
    geom_segment(aes(x = 1, xend = 1.7, y = y_max + 0.09 * y_range, yend = y_max + 0.09 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 1, xend = 1, y = y_max + 0.07 * y_range, yend = y_max + 0.09 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 1.7, xend = 1.7, y = y_max + 0.07 * y_range, yend = y_max + 0.09 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    annotate("text", x = 1.35, y = y_max + 0.13 * y_range, label = pval_T2D_noslgt2_control, size = 4.5) +
    
    geom_segment(aes(x = 1.7, xend = 2.0, y = y_max + 0.03 * y_range, yend = y_max + 0.03 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 1.7, xend = 1.7, y = y_max + 0.01 * y_range, yend = y_max + 0.03 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 2.0, xend = 2.0, y = y_max + 0.01 * y_range, yend = y_max + 0.03 * y_range), 
                 color = "black", size = 0.5, inherit.aes = FALSE) +
    annotate("text", x = 1.85, y = y_max + 0.07 * y_range, label = pval_T2D_comparison, size = 4.5) +
    
    expand_limits(y = y_max + 0.28 * y_range)
  
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
                        results_list[[9]], results_list[[10]],
                        results_list[[11]], results_list[[12]],
                        ncol = 2)

dev.off()

png('C:/Users/netio/Documents/UofW/Rockies/SGLT2ComparisonGroups_FSOC.png', 
    width =1200, height = 1600)  
gridExtra::grid.arrange(results_list[[1]], results_list[[2]], 
                        results_list[[3]], results_list[[4]], 
                        results_list[[5]], results_list[[6]], 
                        results_list[[7]], results_list[[8]],
                        results_list[[9]], results_list[[10]],
                        results_list[[11]], results_list[[12]],
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
                        results_list[[9]], results_list[[10]],
                        results_list[[11]], results_list[[12]],
                        ncol = 2)

dev.off()

png('C:/Users/netio/Documents/UofW/Rockies/SGLT2ComparisonGroups_FSOC_ttest.png', 
    width =1200, height = 1600)  
gridExtra::grid.arrange(results_list[[1]], results_list[[2]], 
                        results_list[[3]], results_list[[4]], 
                        results_list[[5]], results_list[[6]], 
                        results_list[[7]], results_list[[8]],
                        results_list[[9]], results_list[[10]],
                        results_list[[11]], results_list[[12]],
                        ncol = 2)

dev.off()













