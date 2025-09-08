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
  geom_point()+geom_smooth(method='lm', se=F)+
  #geom_smooth(method = 'lm', aes(color = NULL), color = 'black', size = 1.2)+
  theme_classic()+labs(x='aPT Percentage of PT Cells', y = 'dTAL Percentage of TAL Cells')


model_combined <- lm(dTAL_percentage ~ aPT_percentage * group_labels, data = strange_cells)
summary(model_combined)
anova(model_combined)


strange_cells$lean_female_vs_others <- ifelse(strange_cells$group_labels == "Lean_Control_Female", 
                                              "Lean_Control_Female", "Others")

# Test this specific contrast
model_contrast <- lm(dTAL_percentage ~ aPT_percentage * lean_female_vs_others, data = strange_cells)
summary(model_contrast)
anova(model_contrast)

strange_cells$lean_female_vs_others <- ifelse(strange_cells$group_labels == "Lean_Control_Female", 
                                              "Lean_Control_Female", "Others")

# Test this specific contrast
model_contrast <- lm(dTAL_percentage ~ aPT_percentage * lean_female_vs_others, data = strange_cells)
summary(model_contrast)




strange_cells %>% 
  filter(group_labels != 'Lean_Control_Female') %>% 
  with(cor.test(aPT_percentage, dTAL_percentage))
#this is significant! 0.04 p, 0.316 correlation



strange_cells %>% 
  with(cor.test(aPT_percentage, dTAL_percentage))
#this is not significant. p = 0.22, cor = 0.18






#Find relationship between ratios/percentages and clinical variables 


strange_cells <- data.table::fread('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/aPT_dTAL_cellproportions.txt')




harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>%
  arrange(screen_date) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))



dat2 <- dat %>% dplyr::select(record_id, mrn, visit, starts_with('eGFR'), starts_with('fsoc_'),
                              'lc_k2', 'rc_k2', 'lm_k2', 'rm_k2',
                                   'lc_f', 'rc_f', 'lm_f', 'rm_f') %>% filter(visit == 'baseline')


tmp_df <- dat2 %>% dplyr::select(lc_k2, rc_k2, lm_k2, rm_k2,
                                 lc_f, rc_f, lm_f, rm_f)
avg_c_k2 <- tmp_df %>%
  dplyr::select(lc_k2, rc_k2) %>% rowMeans(na.rm=T)

avg_m_k2 <- tmp_df %>% 
  dplyr::select(lm_k2, rm_k2) %>% rowMeans(na.rm=T)

avg_c_f <- tmp_df %>% 
  dplyr::select(lc_f, rc_f) %>% rowMeans(na.rm=T)

avg_m_f <- tmp_df %>% 
  dplyr::select(lm_f, rm_f) %>% rowMeans(na.rm=T)

avg_c_k2_f <- avg_c_k2 / avg_c_f

avg_m_k2_f <- avg_m_k2 / avg_m_f

results <- bind_cols(avg_c_k2, avg_m_k2, avg_c_f, avg_m_f, 
                     avg_c_k2_f, avg_m_k2_f) %>% as.data.frame()
names(results) <- c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 
                    'avg_c_k2_f', 'avg_m_k2_f')

dat2 <- dat2 %>% bind_cols(results)


strange_cells_full <- strange_cells %>% left_join(dat2)

results_df <- data.frame()

for(i in c(11:ncol(strange_cells_full))){
  tmp_var <- names(strange_cells_full)[i]
  tmp_df <- strange_cells_full %>% 
    dplyr::select(group_labels, aPT_percentage, dTAL_percentage, all_of(tmp_var))
  
  names(tmp_df)[4] <- 'Variable'
  tmp_df2 <- tmp_df %>% filter(!is.na(Variable))
  
  # Fix the logical conditions
  n_groups <- length(unique(tmp_df2$group_labels))
  
  if(n_groups < 2){
    next
  } else if(n_groups == 2){
    # Handle 2-group case if needed
    # You can add analysis here for 2-group scenarios
    tmp_df2$group_labels <- tmp_df2$group_labels %>% 
      str_replace(pattern = '_Male', ':Male') %>% 
      str_replace(pattern = '_Female', ':Female')
    tmp_df2 <- tmp_df2 %>% separate(group_labels, into = c('group', 'sex'), sep=':')
    
    for(predictor in c('aPT_percentage', 'dTAL_percentage')){
      if(predictor == 'aPT_percentage'){
        model_full <- lm(Variable ~ aPT_percentage*sex, data = tmp_df2)
        predictor_name <- 'aPT'
      }else{
        model_full <- lm(Variable ~ dTAL_percentage * sex, data = tmp_df2)
        predictor_name <- 'dTAL'
      }
      
      model_summary <- summary(model_full)
      model_anova <- anova(model_full)
      
      temp_results <- data.frame(
        Variable_Name = tmp_var,
        Predictor = predictor_name,
        N_observations = nobs(model_full),
        N_groups = n_groups,
        
        # Model fit
        R_squared = model_summary$r.squared,
        Adj_R_squared = model_summary$adj.r.squared,
        Overall_Model_P = pf(model_summary$fstatistic[1], 
                             model_summary$fstatistic[2], 
                             model_summary$fstatistic[3], 
                             lower.tail = FALSE),
        
        # Main effects (only predictor and sex for 2-group case)
        Predictor_Main_Effect_P = model_anova$`Pr(>F)`[1],      # aPT_percentage or dTAL_percentage
        Group_Main_Effect_P = NA,                               # No group variable
        Sex_Main_Effect_P = model_anova$`Pr(>F)`[2],            # sex
        
        # Two-way interactions (only predictor by sex)
        Predictor_by_Group_P = NA,                              # No group variable
        Predictor_by_Sex_P = model_anova$`Pr(>F)`[3],           # predictor:sex
        Group_by_Sex_P = NA,                                    # No group variable
        
        # Three-way interaction (not possible with only 2 variables)
        Predictor_by_Group_by_Sex_P = NA,
        
        # F-statistics for effect sizes
        Predictor_Main_F = model_anova$`F value`[1],
        Group_Main_F = NA,
        Sex_Main_F = model_anova$`F value`[2],
        Predictor_by_Group_F = NA,
        Predictor_by_Sex_F = model_anova$`F value`[3],
        Group_by_Sex_F = NA,
        Predictor_by_Group_by_Sex_F = NA
      )
      
      results_df <- rbind(results_df, temp_results)
      
    }
    
    
    next
  } else if(n_groups == 4){
    # Parse group_labels into group and sex
    tmp_df2$group_labels <- tmp_df2$group_labels %>% 
      str_replace(pattern = '_Male', ':Male') %>% 
      str_replace(pattern = '_Female', ':Female')
    tmp_df2 <- tmp_df2 %>% separate(group_labels, into = c('group', 'sex'), sep=':')
    tmp_df2$group_labels <- paste0(tmp_df2$group, '_', tmp_df2$sex)
    
    # Loop through both aPT_percentage and dTAL_percentage
    for(predictor in c("aPT_percentage", "dTAL_percentage")){
      
      # Create the appropriate formula
      if(predictor == "aPT_percentage"){
        model_full <- lm(Variable ~ aPT_percentage * group * sex, data = tmp_df2)
        predictor_name <- "aPT"
      } else {
        model_full <- lm(Variable ~ dTAL_percentage * group * sex, data = tmp_df2)
        predictor_name <- "dTAL"
      }
      
      # Fix the summary and anova calls
      model_summary <- summary(model_full)
      model_anova <- anova(model_full)
      
      temp_results <- data.frame(
        Variable_Name = tmp_var,
        Predictor = predictor_name,
        N_observations = nobs(model_full),
        N_groups = n_groups,
        
        # Model fit
        R_squared = model_summary$r.squared,
        Adj_R_squared = model_summary$adj.r.squared,
        Overall_Model_P = pf(model_summary$fstatistic[1], 
                             model_summary$fstatistic[2], 
                             model_summary$fstatistic[3], 
                             lower.tail = FALSE),
        
        # Main effects
        Predictor_Main_Effect_P = model_anova$`Pr(>F)`[1],      # aPT_percentage or dTAL_percentage
        Group_Main_Effect_P = model_anova$`Pr(>F)`[2],         # group  
        Sex_Main_Effect_P = model_anova$`Pr(>F)`[3],           # sex
        
        # Two-way interactions
        Predictor_by_Group_P = model_anova$`Pr(>F)`[4],        # predictor:group
        Predictor_by_Sex_P = model_anova$`Pr(>F)`[5],          # predictor:sex
        Group_by_Sex_P = model_anova$`Pr(>F)`[6],              # group:sex
        
        # Three-way interaction
        Predictor_by_Group_by_Sex_P = model_anova$`Pr(>F)`[7], # predictor:group:sex
        
        # F-statistics for effect sizes
        Predictor_Main_F = model_anova$`F value`[1],
        Group_Main_F = model_anova$`F value`[2],
        Sex_Main_F = model_anova$`F value`[3],
        Predictor_by_Group_F = model_anova$`F value`[4],
        Predictor_by_Sex_F = model_anova$`F value`[5],
        Group_by_Sex_F = model_anova$`F value`[6],
        Predictor_by_Group_by_Sex_F = model_anova$`F value`[7]
      )
      
      results_df <- rbind(results_df, temp_results)
    }
    graph1 <- ggplot(tmp_df, aes(x=group_labels, y=Variable, color = group_labels))+
      geom_boxplot()+geom_point()+theme_classic()+labs(x='Condition Group', y = tmp_var)
    
    graph2 <- ggplot(tmp_df, aes(x=aPT_percentage, y=Variable, color = group_labels))+
      geom_point()+geom_smooth(method='lm', se=F)+theme_classic()+labs(x='aPT Percentage of PT Cells', y = tmp_var)
    
    graph3 <- ggplot(tmp_df, aes(x=dTAL_percentage, y=Variable, color = group_labels))+
      geom_point()+geom_smooth(method='lm', se=F)+theme_classic()+labs(x='dTAL Percentage of TAL Cells', y = tmp_var)
    
    
    png(paste0('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/CellType_Plots/', tmp_var, '_groupdifferences.png'), 
        width= 600)
    print(graph1)
    dev.off()
    
    png(paste0('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/CellType_Plots/', tmp_var, '_aPTdifferences.png'))
    print(graph2)
    dev.off()
    
    png(paste0('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/CellType_Plots/', tmp_var, '_dTALdifferences.png'))
    print(graph3)
    dev.off()
    
  } else {
    next
  }
  
}





write.table(results_df, 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/Threeway_Interactions/FSOC_PET_eGFRResults_aPT_dTAL.txt',
            row.names=F, quote=F, sep='\t')





##Panel analysis: 

plot_data <- results_df %>%
  select(Variable_Name, Predictor, N_groups,
         Predictor_Main_Effect_P, Sex_Main_Effect_P, Group_Main_Effect_P,
         Predictor_by_Sex_P, Predictor_by_Group_P, Group_by_Sex_P,
         Predictor_by_Group_by_Sex_P) %>%
  pivot_longer(cols = ends_with("_P"), 
               names_to = "Effect_Type", 
               values_to = "P_value") %>%
  filter(!is.na(P_value)) %>%  # Remove NA values
  mutate(
    neg_log10_p = -log10(P_value),
    significant = P_value < 0.05,
    Effect_Type = case_when(
      Effect_Type == "Predictor_Main_Effect_P" ~ "Predictor Main Effect",
      Effect_Type == "Sex_Main_Effect_P" ~ "Sex Main Effect", 
      Effect_Type == "Group_Main_Effect_P" ~ "Group Main Effect",
      Effect_Type == "Predictor_by_Sex_P" ~ "Predictor × Sex",
      Effect_Type == "Predictor_by_Group_P" ~ "Predictor × Group",
      Effect_Type == "Group_by_Sex_P" ~ "Group × Sex",
      Effect_Type == "Predictor_by_Group_by_Sex_P" ~ "3-Way Interaction"
    )
  )

# Create multi-panel plot
p1 <- ggplot(plot_data, aes(x = reorder(Variable_Name, neg_log10_p), 
                            y = neg_log10_p, 
                            color = Predictor,
                            shape = as.factor(N_groups))) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7) +
  facet_wrap(~ Effect_Type, scales = "free_x", ncol = 2) +
  coord_flip() +
  scale_color_manual(values = c("aPT" = "#E31A1C", "dTAL" = "#1F78B4")) +
  scale_shape_manual(name = "N Groups", values = c("2" = 16, "4" = 17)) +
  labs(
    title = "Statistical Significance Across All Variables and Effects",
    subtitle = "Red dashed line indicates p = 0.05 threshold",
    x = "Variables",
    y = "-log10(p-value)",
    color = "Predictor Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p1)



#Heatmap 

# Create heatmap of all effects
library(viridis)

heatmap_data <- results_df %>%
  select(Variable_Name, Predictor, 
         Predictor_Main_Effect_P, Sex_Main_Effect_P, Group_Main_Effect_P,
         Predictor_by_Sex_P, Predictor_by_Group_P, Group_by_Sex_P,
         Predictor_by_Group_by_Sex_P) %>%
  pivot_longer(cols = ends_with("_P"), 
               names_to = "Effect_Type", 
               values_to = "P_value") %>%
  filter(!is.na(P_value)) %>%
  mutate(
    neg_log10_p = -log10(P_value),
    Variable_Predictor = paste(Variable_Name, Predictor, sep = "_"),
    Effect_Type = case_when(
      Effect_Type == "Predictor_Main_Effect_P" ~ "Predictor\nMain",
      Effect_Type == "Sex_Main_Effect_P" ~ "Sex\nMain", 
      Effect_Type == "Group_Main_Effect_P" ~ "Group\nMain",
      Effect_Type == "Predictor_by_Sex_P" ~ "Predictor\n× Sex",
      Effect_Type == "Predictor_by_Group_P" ~ "Predictor\n× Group",
      Effect_Type == "Group_by_Sex_P" ~ "Group\n× Sex",
      Effect_Type == "Predictor_by_Group_by_Sex_P" ~ "3-Way\nInteraction"
    )
  )

p3 <- ggplot(heatmap_data, aes(x = Effect_Type, 
                               y = reorder(Variable_Predictor, neg_log10_p),
                               fill = neg_log10_p)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_c(name = "-log10(p)", 
                       option = "plasma",
                       trans = "sqrt") +
  geom_text(aes(label = ifelse(neg_log10_p > -log10(0.05), "*", "")), 
            color = "white", size = 3, fontface = "bold") +
  labs(
    title = "Statistical Significance Heatmap",
    subtitle = "Asterisks (*) indicate p < 0.05",
    x = "Effect Type",
    y = "Variable_Predictor"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(p3)







