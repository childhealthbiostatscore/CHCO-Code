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

pdf("C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/CelltypeDistribution_bySex.pdf")
ggplot(meta.data %>% filter(celltype_rpca != 'ATL'), aes(x=group_labels, fill=celltype_rpca))+
  geom_bar(position='fill')+theme_classic()+labs(x='Condition Group', fill='Cell Type')

dev.off()


so <- subset(so, group %in% c('Type_2_Diabetes', 'Lean_Control'))



















