###Sex Specific Analyses: Investigating aPT/dTAL


#Load in packages for functions we need. 

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







load('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/Line265.RData')


so <- subset(so, group %in% c('Type_2_Diabetes', 'Lean_Control'))

meta.data <- so@meta.data



#aPT analysis
cell_distribution <- meta.data %>% 
  mutate(group_labels = paste0(group, '_', sex)) %>%
  filter(KPMP_celltype %in% c('aPT', 'PT-S1/S2', 'PT-S3')) %>% 
  group_by(record_id, group_labels) %>%
  summarize(
    total_cells = n(), 
    aPT_count = sum(KPMP_celltype == 'aPT'),
    aPT_percentage = (sum(KPMP_celltype == 'aPT') / n()) * 100,
    .groups = 'drop'
  )


pdf("C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/aPT_Percentage_byParticipant.pdf")

ggplot(cell_distribution, aes(x=group_labels, color=group_labels, y= aPT_percentage))+
  geom_point(position=position_jitter(width = 0.1, height = 0))+geom_boxplot()+theme_classic()+labs(x='Condition Group', y='Percent aPT of PT Cells', title = 'aPT Percentage in Each Participant by Group')
dev.off()


pdf("C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/aPT_CountbyParticipant.pdf")

ggplot(cell_distribution, aes(x=group_labels, color=group_labels, y= aPT_count))+
  geom_point(position=position_jitter(width = 0.1, height = 0))+geom_boxplot()+theme_classic()+labs(x='Condition Group', y='Number of aPT Cells', title = 'aPT Cells in Each Participant by Group')
dev.off()


aPT_data <- cell_distribution



#TAL analysis

cell_distribution <- meta.data %>% 
  mutate(group_labels = paste0(group, '_', sex)) %>%
  filter(KPMP_celltype %in% c('dTAL', 'aTAL', 'C-TAL-1', 'C-TAL-2')) %>% 
  group_by(record_id, group_labels) %>%
  summarize(
    total_cells = n(), 
    dTAL_count = sum(KPMP_celltype == 'dTAL'),
    dTAL_percentage = (sum(KPMP_celltype == 'dTAL') / n()) * 100,
    .groups = 'drop'
  )


pdf("C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/dTAL_Percentage_byParticipant.pdf")

ggplot(cell_distribution, aes(x=group_labels, color=group_labels, y= dTAL_percentage))+
  geom_point(position=position_jitter(width = 0.1, height = 0))+geom_boxplot()+theme_classic()+labs(x='Condition Group', y='Percent dTAL in TAL Cells', title = 'dTAL Percentage in Each Participant by Group')
dev.off()


pdf("C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/dTAL_CountbyParticipant.pdf")

ggplot(cell_distribution, aes(x=group_labels, color=group_labels, y= dTAL_count))+
  geom_point(position=position_jitter(width = 0.1, height = 0))+geom_boxplot()+theme_classic()+labs(x='Condition Group', y='Number of dTAL Cells', title = 'dTAL Cells in Each Participant by Group')
dev.off()


TAL_data <- cell_distribution



#Combined for analyses

names(aPT_data) <- c('record_id', 'group_labels', 'PT_totalcells', 'aPT_count', 'aPT_percentage')
names(TAL_data) <- c('record_id', 'group_labels', 'TAL_totalcells', 'dTAL_count', 'dTAL_percentage')

strange_cells <- aPT_data %>% left_join(TAL_data)

write.table(strange_cells, 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/aPT_dTAL_cellproportions.txt', row.names=F, quote=F, sep='\t')

ggplot(strange_cells, aes(x=aPT_percentage, y = dTAL_percentage, color = group_labels))+
  geom_point()+geom_smooth(method='lm', se=F)+geom_smooth(method = 'lm', aes(color = NULL), color = 'black', size = 1.2)+
  theme_classic()+labs(x='aPT Percentage of PT Cells', y = 'dTAL Percentage of TAL Cells')





strange_cells %>% 
  filter(group_labels != 'Lean_Control_Female') %>% 
  with(cor.test(aPT_percentage, dTAL_percentage))
#this is significant! 0.04 p, 0.316 correlation



strange_cells %>% 
  with(cor.test(aPT_percentage, dTAL_percentage))
#this is not significant. p = 0.22, cor = 0.18




