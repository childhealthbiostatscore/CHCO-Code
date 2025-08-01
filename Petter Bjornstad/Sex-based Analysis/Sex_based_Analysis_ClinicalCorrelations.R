#Sex-Specific Analyses of Clinical Correlates

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



harmonized_data <- read.csv("C:/Users/netio/Documents/Harmonized_data/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(screen_date) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


dat <- dat %>% filter(group %in% c('Lean Control', 'Type 2 Diabetes'))


library(ggplot2)

names.dat <- names(dat)


ggplot(dat %>% filter(acr_u < 2000 & diabetes_duration < 15) %>% 
         filter(!is.na(sex)), 
       aes(x=diabetes_duration, y = acr_u, color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()



ggplot(dat %>% filter(acr_u < 2000 & diabetes_duration < 15) %>% 
         filter(!is.na(sex)), 
       aes(x=diabetes_duration, y = log10(acr_u), color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()

ggplot(dat %>% filter(acr_u < 2000 & diabetes_duration < 15) %>% 
         filter(!is.na(sex)), 
       aes(x=diabetes_duration, y = log10(acr_u), color = sex))+
  geom_density_2d()+
  theme_classic()





ggplot(dat %>% filter(acr_u < 2000 & diabetes_duration < 15) %>% 
         filter(!is.na(sex)), 
       aes(x=diabetes_duration, y = log10(acr_u), color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()




gfr_variables <- names.dat[str_which(names.dat, 'gfr')]

pdf('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/Clinical/GFR_graphs.pdf')
for(i in c(1:length(gfr_variables))){
  index_var <- which(names(dat) == gfr_variables[i])
  tmp <- dat %>% select(diabetes_duration, sex, index_var)
  names(tmp)[3] <- 'Variable'
  
  gfr_graphs <- ggplot(tmp %>% filter(diabetes_duration < 15) %>% 
           filter(!is.na(sex)), 
         aes(x=diabetes_duration, y = Variable, color = sex))+
    geom_point()+
    geom_smooth(method = 'lm')+
    theme_classic()+labs(y=gfr_variables[i])
  print(gfr_graphs)
  
  
}
dev.off()













