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





tests <- c('fsoc_l_cortex', 'fsoc_r_cortex', 
           'fsoc_l_kidney', 'fsoc_r_kidney', 
           'fsoc_l_medulla', 'fsoc_r_medulla')


graphs <- list()




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
                fsoc_l_kidney, fsoc_r_kidney, fsoc_l_medulla, fsoc_r_medulla)


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



pdf('C:/Users/netio/Documents/UofW/Rockies/SGLT2ComparisonGroups_FSOC.pdf', 
    width =20, height = 20)  
gridExtra::grid.arrange(results_list[[1]], results_list[[2]], 
                        results_list[[3]], results_list[[4]], 
                        results_list[[5]], results_list[[6]], ncol = 2)

dev.off()













