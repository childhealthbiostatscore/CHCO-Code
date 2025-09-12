### Module Score Analysis
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


#Lean Control vs. T2D Score


load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')

so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

so_subset <- AddModuleScore(object = so_subset, 
                            features = list(tca_genes), 
                            name = 'TCA_score')

so_subset <- AddModuleScore(object = so_subset, 
                            features = list(ox_phos_genes),
                            name = 'OxPhos_score')

meta.data <- so_subset@meta.data




VlnPlot(so_subset, features = 'TCA_score1', group.by = 'group',
        pt.size = 0)

library(ggplot2)

# Extract the data
plot_data <- data.frame(
  TCA_score = so_subset$TCA_score1,
  OxPhos_Score = so_subset$OxPhos_score1,
  condition = so_subset$group
)

ggplot(plot_data, aes(x = condition, y = TCA_score, fill = condition)) +
  geom_violin(alpha = 0.7) +                    # Violin plot
  geom_boxplot(width = 0.2, alpha = 0.8) + 
  scale_x_discrete(labels = c("Lean_Control" = "Lean Control", 
                              "Type_2_Diabetes" = "Type 2 Diabetes (no SGLT2)")) +
  scale_fill_discrete(labels = c("Lean_Control" = "Lean Control", 
                                 "Type_2_Diabetes" = "Type 2 Diabetes (no SGLT2)")) +
  # Narrow boxplot on top
  theme_classic() +
  labs(y = "TCA Module Score", x = "Condition") # Match your original colors


write.table(plot_data, 'C:/Users/netio/Documents/UofW/Rockies/Module_scores/LC_vs_T2D_noSGLT2_TCA_OxPhos_ModuleScores.txt', 
            row.names=F, quote=F, sep='\t')

#T2D Score 
remove(list=ls())



load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')


load("C:/Users/netio/Downloads/TCA_genes.txt")
load('C:/Users/netio/Downloads/OxPhos_genes.txt')

so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

test <- so_subset@meta.data

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


dat <- dat %>% dplyr::select(record_id, epic_sglti2_1) %>% 
  filter(!is.na(epic_sglti2_1))

test$epic_sglti2_1 <- NULL

test <- test %>% left_join(dat, by='record_id')

so_subset@meta.data$epic_sglti2_1 <- test$epic_sglti2_1





so_subset <- AddModuleScore(object = so_subset, 
                            features = list(tca_genes), 
                            name = 'TCA_score')

so_subset <- AddModuleScore(object = so_subset, 
                            features = list(ox_phos_genes),
                            name = 'OxPhos_score')

meta.data <- so_subset@meta.data




VlnPlot(so_subset, features = 'TCA_score1', group.by = 'group',
        pt.size = 0)

library(ggplot2)

# Extract the data
plot_data <- data.frame(
  TCA_score = so_subset$TCA_score1,
  OxPhos_Score = so_subset$OxPhos_score1,
  condition = so_subset$group
)

ggplot(plot_data, aes(x = condition, y = TCA_score, fill = condition)) +
  geom_violin(alpha = 0.7) +                    # Violin plot
  geom_boxplot(width = 0.2, alpha = 0.8) + 
  scale_x_discrete(labels = c("Lean_Control" = "Lean Control", 
                              "Type_2_Diabetes" = "Type 2 Diabetes (no SGLT2)")) +
  scale_fill_discrete(labels = c("Lean_Control" = "Lean Control", 
                                 "Type_2_Diabetes" = "Type 2 Diabetes (no SGLT2)")) +
  # Narrow boxplot on top
  theme_classic() +
  labs(y = "TCA Module Score", x = "Condition") # Match your original colors


write.table(plot_data, 'C:/Users/netio/Documents/UofW/Rockies/Module_scores/T2D_SGLT2_noSGLT2_TCA_OxPhos_ModuleScores.txt', 
            row.names=F, quote=F, sep='\t')










































