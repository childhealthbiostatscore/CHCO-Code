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




#VlnPlot(so_subset, features = 'TCA_score1', group.by = 'group',
#        pt.size = 0)

library(ggplot2)

# Extract the data
plot_data <- data.frame(
  record_id = so_subset$record_id,
  mrn = so_subset$mrn, 
  TCA_score = so_subset$TCA_score1,
  OxPhos_Score = so_subset$OxPhos_score1,
  condition = so_subset$group,
  KPMP_celltype = so_subset$KPMP_celltype,
  celltype2 = so_subset$celltype2
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





so_subset$celltype1 <- case_when(grepl("PT-",so_subset$celltype_rpca)~"PT",
                                 grepl("TAL-",so_subset$celltype_rpca)~"TAL",
                                 grepl("EC-",so_subset$celltype_rpca)~"EC",
                                 grepl("POD",so_subset$celltype_rpca)~"POD",
                                 grepl("MAC",so_subset$celltype_rpca)~"MAC",
                                 grepl("MON",so_subset$celltype_rpca)~"MON",
                                 grepl("PC-",so_subset$celltype_rpca)~"PC",
                                 grepl("FIB",so_subset$celltype_rpca)~"FIB_MC_VSMC",
                                 grepl("DTL",so_subset$celltype_rpca)~"DTL",
                                 so_subset$celltype_rpca=="DCT"~"DCT",
                                 so_subset$celltype_rpca=="ATL"~"ATL",
                                 so_subset$celltype_rpca=="B"~"B",
                                 so_subset$celltype_rpca=="T"~"T")
so_subset$celltype1 <- as.character(so_subset$celltype1)

so_subset$KPMP_celltype2 <- as.character(so_subset$KPMP_celltype)
so_subset$celltype2 <- ifelse(so_subset$KPMP_celltype=="aPT" | 
                                so_subset$KPMP_celltype=="PT-S1/S2" | 
                                so_subset$KPMP_celltype == "PT-S3","PT",
                              ifelse(grepl("TAL",so_subset$KPMP_celltype),"TAL",
                                     ifelse(grepl("EC-",so_subset$KPMP_celltype),"EC",so_subset$KPMP_celltype2)))


so_subset$DCT_celltype <- ifelse((so_subset$KPMP_celltype=="DCT" | 
                                    so_subset$KPMP_celltype=="dDCT"), "DCT","Non-DCT")








so_subset <- AddModuleScore(object = so_subset, 
                            features = list(tca_genes), 
                            name = 'TCA_score')

so_subset <- AddModuleScore(object = so_subset, 
                            features = list(ox_phos_genes),
                            name = 'OxPhos_score')

meta.data <- so_subset@meta.data




#VlnPlot(so_subset, features = 'TCA_score1', group.by = 'group',
#        pt.size = 0)

library(ggplot2)

# Extract the data
plot_data <- data.frame(
  record_id = so_subset$record_id,
  mrn = so_subset$mrn, 
  TCA_score = so_subset$TCA_score1,
  OxPhos_Score = so_subset$OxPhos_score1,
  condition = so_subset$group, 
  sglt2 = so_subset$epic_sglti2_1,
  KPMP_celltype = so_subset$KPMP_celltype,
  celltype2 = so_subset$celltype2
)

write.table(plot_data, 'C:/Users/netio/Documents/UofW/Rockies/Module_scores/T2D_SGLT2_noSGLT2_TCA_OxPhos_ModuleScores.txt', 
            row.names=F, quote=F, sep='\t')







## Combinations and plots

remove(list=ls())


lc_t2d <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Module_scores/LC_vs_T2D_noSGLT2_TCA_OxPhos_ModuleScores.txt')
t2d_sglt2 <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Module_scores/T2D_SGLT2_noSGLT2_TCA_OxPhos_ModuleScores.txt')

lc_t2d <- lc_t2d %>% 
  mutate(condition = ifelse(condition == 'Lean_Control', 'Lean Control', 'T2D-No SGLT2'))

t2d_all <- t2d_sglt2 %>% 
  dplyr::select(record_id, mrn, TCA_score, OxPhos_Score, KPMP_celltype, celltype2)
t2d_all$condition <- 'T2D Combined'

t2d_sglt2 <- t2d_sglt2 %>% 
  anti_join(lc_t2d, by='record_id') %>%
  mutate(condition = ifelse(sglt2 == 'Yes', 'T2D-SGLT2', 'T2D-No SGLT2')) %>% 
  dplyr::select(-sglt2)


full_results <- bind_rows(list(lc_t2d, t2d_all, t2d_sglt2))
full_results$group2 <- full_results$condition



#boxplot function 
boxplot_function <- function(data, variable, label, method, celltype){
  if(celltype == 'All'){
    data <- data
  }else if(celltype %in% c('PT', 'TAL', 'EC')){
    data <- data %>% filter(celltype2 == celltype)
  }else{
    data <- data %>% filter(KPMP_celltype2 == celltype)
  }
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  
  
  var_index <- which(names(data) == variable)
  data <- data %>% dplyr::select(group2, var_index)
  names(data)[2] <- 'Variable'
  
  data <- data %>% mutate(position = ifelse(group2 == 'T2D-No SGLT2', 1.7, ifelse(group2 == 'T2D-SGLT2', 2.0, NA)))
  
  
  plot <- ggplot(data %>% dplyr::filter(group2 %in% c('Lean Control', 'T2D Combined')), aes(x = group2, y = Variable, fill = group2))  +
    geom_boxplot(width = 1.3, size = 1)+
    scale_fill_manual(values = c("#c2dfe3", "#fff9ec", "#fcb1a6", "#fb6376")) +
    geom_boxplot(data = data %>%
                   dplyr::filter(group2 %in% c('T2D-No SGLT2', 'T2D-SGLT2')), 
                 aes(x = position, y=Variable, fill = group2), width = 0.1, size = 1)+
    labs(x= 'Study Group', y = label, fill = 'Study Group', title = paste0('Analysis of ', label, ' in ', 
                                                                           celltype2, ' Cells'))+
    theme_minimal()+
    theme(axis.text.x = element_blank(),
          text = element_text(size = 20))
  
  if(method == 'ANOVA'){
 #   model <- aov(Variable ~ group2, data = data)
#    model_results <- TukeyHSD(model, conf.level = 0.95)$group2 %>% 
#      as.data.frame()
#    
#    model_results <- model_results %>% 
#      mutate(pvalue = ifelse(`p adj` < 0.001, '< 0.001', 
#                             paste0('p = ', round(`p adj`, 3))))
#    
#    pval_T2D_noslgt2_control <- model_results$pvalue[which(rownames(model_results) == 'T2D-No SGLT2-Lean Control')] %>%
#      as.character()
#    pval_T2D_total_control <- model_results$pvalue[which(rownames(model_results) == 'T2D Combined-Lean Control')] %>% 
#      as.character()
#    pval_T2D_comparison <- model_results$pvalue[which(rownames(model_results) == 'T2D-SGLT2-T2D-No SGLT2')] %>% 
#      as.character()
  }else if(method == 't-test'){
    
    
 #   tmp <- data %>% filter(group2 %in% c('T2D-No SGLT2', 'Lean Control'))
#    model1 <- t.test(Variable ~ group2, data = tmp)
#    pval_T2D_noslgt2_control <- ifelse(model1$p.value < 0.001, '< 0.001',
#                                       paste0('p = ', round(model1$p.value, 3)))
#    
#    tmp <- data %>% filter(group2 %in% c('T2D Combined', 'Lean Control'))
#    model1 <- t.test(Variable ~ group2, data = tmp)
#    pval_T2D_total_control <- ifelse(model1$p.value < 0.001, '< 0.001',
#                                     paste0('p = ', round(model1$p.value, 3)))
#    
#    tmp <- data %>% filter(group2 %in% c('T2D-No SGLT2', 'T2D-SGLT2'))
#    model1 <- t.test(Variable ~ group2, data = tmp)
#    pval_T2D_comparison <- ifelse(model1$p.value < 0.001, '< 0.001',
#                                  paste0('p = ', round(model1$p.value, 3)))
    
    
  }
  
#  y_max <- max(data$Variable, na.rm = TRUE)
#  y_range <- diff(range(data$Variable, na.rm = TRUE))
#  
#  
#  plot <- plot + 
#    geom_segment(aes(x = 1, xend = 2, y = y_max + 0.15 * y_range, yend = y_max + 0.15 * y_range), 
#                 color = "black", size = 0.5, inherit.aes = FALSE) +
#    geom_segment(aes(x = 1, xend = 1, y = y_max + 0.13 * y_range, yend = y_max + 0.15 * y_range), 
#                 color = "black", size = 0.5, inherit.aes = FALSE) +
#    geom_segment(aes(x = 2, xend = 2, y = y_max + 0.13 * y_range, yend = y_max + 0.15 * y_range), 
#                 color = "black", size = 0.5, inherit.aes = FALSE) +
#    annotate("text", x = 1.5, y = y_max + 0.19 * y_range, label = pval_T2D_total_control, size = 4.5) +
#    
#    geom_segment(aes(x = 1, xend = 1.7, y = y_max + 0.09 * y_range, yend = y_max + 0.09 * y_range), 
##                 color = "black", size = 0.5, inherit.aes = FALSE) +
#    geom_segment(aes(x = 1, xend = 1, y = y_max + 0.07 * y_range, yend = y_max + 0.09 * y_range), 
#                 color = "black", size = 0.5, inherit.aes = FALSE) +
#    geom_segment(aes(x = 1.7, xend = 1.7, y = y_max + 0.07 * y_range, yend = y_max + 0.09 * y_range), 
#                 color = "black", size = 0.5, inherit.aes = FALSE) +
#    annotate("text", x = 1.35, y = y_max + 0.13 * y_range, label = pval_T2D_noslgt2_control, size = 4.5) +
#    
#    geom_segment(aes(x = 1.7, xend = 2.0, y = y_max + 0.03 * y_range, yend = y_max + 0.03 * y_range), 
#                 color = "black", size = 0.5, inherit.aes = FALSE) +
#    geom_segment(aes(x = 1.7, xend = 1.7, y = y_max + 0.01 * y_range, yend = y_max + 0.03 * y_range), 
#                 color = "black", size = 0.5, inherit.aes = FALSE) +
#    geom_segment(aes(x = 2.0, xend = 2.0, y = y_max + 0.01 * y_range, yend = y_max + 0.03 * y_range), 
#                 color = "black", size = 0.5, inherit.aes = FALSE) +
#    annotate("text", x = 1.85, y = y_max + 0.07 * y_range, label = pval_T2D_comparison, size = 4.5) +
#    
#    expand_limits(y = y_max + 0.28 * y_range)
  
  print(plot)
  
  
  
  
}


graph_all <- boxplot_function(data = full_results, variable = 'TCA_score', label = 'TCA Module Score', method = 'ttest', celltype = 'All')






### Heatmaps
library(dplyr)
library(ggplot2)
library(tidyr)

# Calculate mean scores by cell type and condition
heatmap_data <- full_results %>% 
  filter(KPMP_celltype %in% c('C-TAL-1', 'C-TAL-2', 'dTAL', 'aPT', 'PT-S1/S2', 'PT-S3')) %>% 
  group_by(KPMP_celltype, condition) %>%
  summarise(
    TCA_mean = mean(TCA_score, na.rm = TRUE),
    OxPhos_mean = mean(OxPhos_Score, na.rm = TRUE),
    .groups = 'drop'
  )

print(heatmap_data)

# Step 2: Reshape for Heatmap
# Convert to long format
heatmap_long <- heatmap_data %>%
  pivot_longer(cols = c(TCA_mean, OxPhos_mean), 
               names_to = "score_type", 
               values_to = "mean_score") %>%
  mutate(score_type = gsub("_mean", "", score_type)) %>%
  # Reorder the cell types to group TAL and PT together
  mutate(KPMP_celltype = factor(KPMP_celltype, 
                                levels = c('C-TAL-1','C-TAL-2', 'dTAL', 
                                           'PT-S1/S2', 'PT-S3', 'aPT')))

# Step 3: Create the Heatmap
tmp_plot <- ggplot(heatmap_long, aes(x = condition, y = KPMP_celltype, fill = mean_score)) +
  geom_tile(color = "white", size = 0.5) +
  facet_wrap(~score_type, scales = "free") +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "salmon", 
                       midpoint = 0, name = "Mean Score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(x = "Condition", y = "Cell Type", title = "TCA and OxPhos Module Scores")

pdf('C:/Users/netio/Downloads/PT_TAL_TCA_OxPhos_ModuleScores.pdf', width = 15, height = 15)
tmp_plot
dev.off()




### Heatmaps
library(dplyr)
library(ggplot2)
library(tidyr)

# Calculate mean scores by cell type and condition
heatmap_data <- full_results %>% 
  filter(KPMP_celltype %in% c('C-TAL-1', 'C-TAL-2', 'dTAL', 'aPT', 'PT-S1/S2', 'PT-S3')) %>% 
  group_by(KPMP_celltype, condition) %>%
  summarise(
    TCA_mean = mean(TCA_score, na.rm = TRUE),
    OxPhos_mean = mean(OxPhos_Score, na.rm = TRUE),
    .groups = 'drop'
  )

print(heatmap_data)

# Step 2: Reshape for Heatmap
# Convert to long format
heatmap_long <- heatmap_data %>%
  filter(condition %in% c('Lean Control', 'T2D Combined')) %>% 
  pivot_longer(cols = c(TCA_mean, OxPhos_mean), 
               names_to = "score_type", 
               values_to = "mean_score") %>%
  mutate(score_type = gsub("_mean", "", score_type)) %>%
  # Reorder the cell types to group TAL and PT together
  mutate(KPMP_celltype = factor(KPMP_celltype, 
                                levels = c('C-TAL-1','C-TAL-2', 'dTAL', 
                                           'PT-S1/S2', 'PT-S3', 'aPT')))

# Step 3: Create the Heatmap
tmp_plot <- ggplot(heatmap_long, aes(x = condition, y = KPMP_celltype, fill = mean_score)) +
  geom_tile(color = "white", size = 0.5) +
  facet_wrap(~score_type, scales = "free") +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "salmon", 
                       midpoint = 0, name = "Mean Score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(x = "Condition", y = "Cell Type", title = "TCA and OxPhos Module Scores")

pdf('C:/Users/netio/Downloads/heatmap_LCvsT2D.pdf', width = 15, height = 15)
tmp_plot
dev.off()












