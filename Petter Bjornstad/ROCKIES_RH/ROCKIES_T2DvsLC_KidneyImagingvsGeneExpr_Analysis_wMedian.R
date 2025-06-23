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












fixed_data <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/T2D_rawimagingdata.txt')

# #Calculate K2 and F variables
fixed_data <- fixed_data %>%
  mutate(avg_c_k2 = (lc_k2+rc_k2)/2) %>%
  mutate(avg_m_k2 = (lm_k2+rm_k2)/2) %>%
  mutate(avg_c_f = (lc_f+rc_f)/2) %>%
  mutate(avg_m_f = (lm_f+rm_f)/2)
fixed_data <- fixed_data %>%
  rowwise() %>%
  mutate(avg_c_k2_f = (avg_c_k2/avg_c_f)) %>%
  mutate(avg_m_k2_f = (avg_m_k2/avg_m_f))

med_c_k2 <- median(fixed_data$avg_c_k2, na.rm=T)
med_m_k2 <- median(fixed_data$avg_m_k2, na.rm=T)
med_c_f <- median(fixed_data$avg_c_f, na.rm=T)
med_m_f <- median(fixed_data$avg_m_f, na.rm=T)
med_c_k2_f <- median(fixed_data$avg_c_k2_f, na.rm=T)
med_m_k2_f <- median(fixed_data$avg_m_k2_f, na.rm=T)


fixed_data <- fixed_data %>% 
  rowwise() %>% 
  mutate(avg_c_k2_med = ifelse(avg_c_k2 >= med_c_k2, 'Above Median', 'Below Median'),
         avg_m_k2_med = ifelse(avg_m_k2 >= med_m_k2, 'Above Median', 'Below Median'),
         avg_c_f_med = ifelse(avg_c_f >= med_c_f, 'Above Median', 'Below Median'),
         avg_m_f_med = ifelse(avg_m_f >= med_m_f, 'Above Median', 'Below Median'),
         avg_c_k2_f_med = ifelse(avg_c_k2_f >= med_c_k2_f, 'Above Median', 'Below Median'),
         avg_m_k2_f_med = ifelse(avg_m_k2_f >= med_m_k2_f, 'Above Median', 'Below Median'))



#Try to add this data from dat 
so_subset@meta.data <- so_subset@meta.data[, !colnames(so_subset@meta.data) %in% c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 
                                                                                   'avg_m_f', 
                                                                                   'avg_c_k2_f', 'avg_m_k2_f',
                                                                                   'lc_k2', 'rc_k2', 'lm_k2', 
                                                                                   'rm_k2', 'lm_f', 'rm_f',
                                                                                   'avg_c_k2_med', 'avg_m_k2_med',
                                                                                   'avg_c_f_med', 'avg_m_f_med', 
                                                                                   'avg_c_k2_f_med', 'avg_m_k2_f_med')]

so_subset@meta.data <- so_subset@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  left_join(fixed_data, by = "mrn") %>%
  tibble::column_to_rownames("cell_id")


so_subset$avg_c_k2_med <- factor(so_subset$avg_c_k2_med)
so_subset$avg_c_k2_med <- relevel(so_subset$avg_c_k2_med,"Below Median")

so_subset$avg_m_k2_med <- factor(so_subset$avg_m_k2_med)
so_subset$avg_m_k2_med <- relevel(so_subset$avg_m_k2_med,"Below Median")

so_subset$avg_c_f_med <- factor(so_subset$avg_c_f_med)
so_subset$avg_c_f_med <- relevel(so_subset$avg_c_f_med,"Below Median")

so_subset$avg_m_f_med <- factor(so_subset$avg_m_f_med)
so_subset$avg_m_f_med <- relevel(so_subset$avg_m_f_med,"Below Median")

so_subset$avg_c_k2_f_med <- factor(so_subset$avg_c_k2_f_med)
so_subset$avg_c_k2_f_med <- relevel(so_subset$avg_c_k2_f_med,"Below Median")

so_subset$avg_m_k2_f_med <- factor(so_subset$avg_m_k2_f_med)
so_subset$avg_m_k2_f_med <- relevel(so_subset$avg_m_k2_f_med,"Below Median")


















