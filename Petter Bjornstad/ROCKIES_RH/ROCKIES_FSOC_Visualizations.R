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
  filter(study %in% c('RENAL-HEIR', 'RENAL-HEIRitage') | record_id == 'IT_19') %>%
  filter(group != 'Obese Control') %>% 
  dplyr::select(mrn, record_id, visit, group, group2, starts_with('fsoc'), bmi, 
                epic_sglti2_1, age, sex, epic_mfm_1, epic_insulin_1, epic_glp1ra_1)

table1::table1(~age + sex + bmi+epic_mfm_1+epic_insulin_1+epic_glp1ra_1+epic_sglti2_1 | group2,data=dat2 %>% 
                 filter(!is.na(group2)))


tests <- c('fsoc_l_cortex', 'fsoc_r_cortex', 
           'fsoc_l_kidney', 'fsoc_r_kidney', 
           'fsoc_l_medulla', 'fsoc_r_medulla', 
           'fsoc_l_combined', 'fsoc_r_combined', 
           'fsoc_full_combined')


graphs <- list()


dat_results <- dat2 %>% mutate(fsoc_l_combined = (fsoc_l_cortex +fsoc_l_kidney+ fsoc_l_medulla)/3,
                               fsoc_r_combined = (fsoc_r_cortex +fsoc_r_kidney +fsoc_r_medulla)/3,
                               fsoc_full_combined = (fsoc_l_cortex + fsoc_r_cortex +
                                                         fsoc_l_kidney + fsoc_r_kidney + 
                                                         fsoc_l_medulla + fsoc_r_medulla)/6)



medications <- readxl::read_xlsx("C:/Users/netio/Documents/Harmonized_data/Biopsies_w_mrn_Oct3.xlsx")
medications <- medications %>% dplyr::select(mrn, ends_with('_1'), -starts_with('ever_'))
names(medications) <- str_replace(names(medications), pattern = '_1', replacement = '')



need_med_info <- dat_results %>% filter(group == 'Type 2 Diabetes') %>% filter(is.na(group2))

med_small <- medications %>% filter(mrn %in% need_med_info$mrn)

for(i in c(1:nrow(med_small))){
  if(med_small$sglti2[i] == 'No'){
    dat_results$group2[which(dat_results$mrn == med_small$mrn[i])] <- 'T2D-No SGLTi2'
  }else if(med_small$sglti2[i] == 'Yes'){
    dat_results$group2[which(dat_results$mrn == med_small$mrn[i])] <- 'T2D-SGLTi2'    
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
    dat_results$group2[which(dat_results$record_id == RH_small$Subject[i])] <- 'T2D-No SGLTi2'
  }else if(RH_small$SGLT2[i] == 'Yes'){
    dat_results$group2[which(dat_results$record_id == RH_small$Subject[i])] <- 'T2D-SGLTi2'    
  }else{
    next
  }
  
}

for(i in c(1:nrow(RH2_small))){
  if(RH2_small$SGLT2[i] == 'No'){
    dat_results$group2[which(dat_results$mrn == RH2_small$mrn[i])] <- 'T2D-No SGLTi2'
  }else if(RH2_small$SGLT2[i] == 'Yes'){
    dat_results$group2[which(dat_results$mrn == RH2_small$mrn[i])] <- 'T2D-SGLTi2'    
  }else{
    next
  }
  
}





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
  mutate(group2 = factor(group2, levels = c('Lean Control', 
                                            'T2D-No SGLTi2', 
                                            'T2D-SGLTi2',
                                            'T2D Combined')))


df_plot <- df_plot %>% 
  dplyr::select(mrn, group, group2, fsoc_l_cortex, fsoc_r_cortex, 
                fsoc_l_kidney, fsoc_r_kidney, fsoc_l_medulla, fsoc_r_medulla,
                fsoc_l_combined, fsoc_r_combined, fsoc_full_combined)


df_plot <- df_plot %>% mutate(position = ifelse(group2 == 'T2D-No SGLTi2', 1.7, ifelse(group2 == 'T2D-SGLTi2', 2.0, NA)))



boxplot_function <- function(data, variable, label){
  
  var_index <- which(names(data) == variable)
  data <- data %>% dplyr::select(group2, var_index)
  names(data)[2] <- 'Variable'
  
  data <- data %>% mutate(position = ifelse(group2 == 'T2D-No SGLTi2', 1.7, ifelse(group2 == 'T2D-SGLTi2', 2.0, NA)))
  
  
  ggplot(data %>% dplyr::filter(group2 %in% c('Lean Control', 'T2D Combined')), aes(x = group2, y = Variable, fill = group2))  +
    geom_boxplot(width = 1.3, size = 1)+
    scale_fill_manual(values = c("#c2dfe3", "#fff9ec", "#fcb1a6", "#fb6376")) +
    geom_boxplot(data = data %>%
                   dplyr::filter(group2 %in% c('T2D-No SGLTi2', 'T2D-SGLTi2')), 
                 aes(x = position, y=Variable, fill = group2), width = 0.1, size = 1)+
    labs(x= 'Study Group', y = label, fill = 'Study Group')+
    theme_minimal()+
    theme(axis.text.x = element_blank(),
          text = element_text(size = 20))
  
}






for(i in c(1:length(tests))){
  if(i == 1){
    results_list <- list()
  }
  results_list[[i]]<- boxplot_function(df_plot, tests[i], tests[i])
}


library(grid)

caption <- grid::textGrob(paste0("Maximum Participants in ", max_num_trait, ": LC: ", max_num_lc, 
                           ';T2D-noSGLT2: ', max_num_nosglt2, ';T2D-SGLT2: ', max_num_sglt2, 
                           'Minimum Participants in Combined Traits: LC: ', min_num_lc, 
                           ';T2D-noSGLT2: ', min_num_nosglt2, ';T2D-SGLT2: ', min_num_sglt2), gp = gpar(fontsize = 12, fontface = "italic"))


pdf('C:/Users/netio/Documents/UofW/Rockies/SGLT2ComparisonGroups_FSOC.pdf', 
    width =20, height = 20)  
gridExtra::grid.arrange(results_list[[1]], results_list[[2]], 
                        results_list[[3]], results_list[[4]], 
                        results_list[[5]], results_list[[6]], 
                        results_list[[7]], results_list[[8]],
                        results_list[[9]], ncol = 2, 
                        bottom = caption)

dev.off()













