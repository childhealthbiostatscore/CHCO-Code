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

harmonized_data <- read.csv("C:/Users/netio/Documents/Harmonized_data/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(screen_date) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


PET_avg <- function(data){
  tmp_df <- data %>% dplyr::select(lc_k2, rc_k2, lm_k2, rm_k2,
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
  
  return(results)
  
}


tmp_results <- PET_avg(dat)


dat_results <- dat %>% dplyr::select(-avg_c_k2, -avg_m_k2, -avg_c_f)

dat_results <- dat_results %>% bind_cols(tmp_results)

dat_results <- dat_results %>% filter(!is.na(avg_c_k2))

table1::table1(~age + sex + bmi + study + group + epic_sglti2_1 + 
                 avg_c_k2 + avg_c_f + avg_m_k2 + avg_m_f + avg_c_k2_f + avg_m_k2_f | group, data = dat_results)










tests <- c('avg_c_k2', 'avg_c_f', 
           'avg_m_k2', 'avg_m_f', 
           'avg_c_k2_f', 'avg_m_k2_f')


graphs <- list()
tests <- c('avg_c_k2', 'avg_c_f', 
           'avg_m_k2', 'avg_m_f', 
           'avg_c_k2_f', 'avg_m_k2_f')



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
  dplyr::select(mrn, group, group2, avg_c_k2, avg_c_f, 
                avg_m_k2, avg_m_f, avg_c_k2_f, avg_m_k2_f)


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
  
  
  
pdf('C:/Users/netio/Documents/UofW/Rockies/SGLT2ComparisonGroups_KidneyImaging.pdf', 
    width =20, height = 20)  
gridExtra::grid.arrange(results_list[[1]], results_list[[2]], 
                        results_list[[3]], results_list[[4]], 
                        results_list[[5]], results_list[[6]], ncol = 2)

dev.off()

  
  
  
  
  










  ggplot(subset(boxplot_adpkd_combined_dat, 
                (pet_param == "avg_f"|pet_param == "avg_k1") & 
                  group_w_class != "ADPKD"),
         mapping = aes(x = group,
                       y = value,
                       fill = group))  +
  scale_fill_manual(values = c("#c2dfe3", "#fff9ec", 
                               "#fcb1a6","#fb6376")) +
  labs(x = NULL,
       y = NULL,
       fill = NULL) +
  geom_boxplot(outlier.shape = NA) +
  geom_rect(data=data.frame(pet_param="avg_f"), 
            aes(ymin = dat_quantile$mild_mod_f[2], 
                ymax = dat_quantile$mild_mod_f[4], 
                xmin = 1.75, 
                xmax = 1.85),
            fill = "#fcb1a6", alpha = 0.9, inherit.aes = F) +
  geom_rect(data=data.frame(pet_param="avg_f"), 
            aes(ymin = dat_quantile$severe_f[2], 
                ymax = dat_quantile$severe_f[4], 
                xmin = 2.15, 
                xmax = 2.25),
            fill = "#fb6376", alpha = 0.9, inherit.aes = F) +
  geom_segment(data=data.frame(pet_param="avg_f"), 
               aes(y = dat_quantile$mild_mod_f[3],
                   yend = dat_quantile$mild_mod_f[3],
                   x = 1.75, xend = 1.85), linewidth = 0.4, inherit.aes = F) +
  geom_segment(data=data.frame(pet_param="avg_f"), 
               aes(y = dat_quantile$severe_f[3],
                   yend = dat_quantile$severe_f[3],
                   x = 2.15, xend = 2.25), linewidth = 0.4, inherit.aes = F) +
  geom_segment(data=data.frame(pet_param="avg_f"), 
               aes(y = dat_quantile$mild_mod_f[1],
                   yend = dat_quantile$mild_mod_f[5],
                   x = 1.8, xend = 1.8),
               linetype = "dashed", linewidth = 0.2, inherit.aes = F) +
  geom_segment(data=data.frame(pet_param="avg_f"), 
               aes(y = dat_quantile$severe_f[1],
                   yend = dat_quantile$severe_f[5],
                   x = 2.2, xend = 2.2),
               linetype = "dashed", linewidth = 0.2, inherit.aes = F) +
  geom_rect(data=data.frame(pet_param="avg_k1"), 
            aes(ymin = dat_quantile$mild_mod_k1[2], 
                ymax = dat_quantile$mild_mod_k1[4], 
                xmin = 1.75, 
                xmax = 1.85),
            fill = "#fcb1a6", alpha = 0.9, inherit.aes = F) +
  geom_rect(data=data.frame(pet_param="avg_k1"), 
            aes(ymin = dat_quantile$severe_k1[2], 
                ymax = dat_quantile$severe_k1[4], 
                xmin = 2.15, 
                xmax = 2.25),
            fill = "#fb6376", alpha = 0.9, inherit.aes = F) +
  geom_segment(data=data.frame(pet_param="avg_k1"), 
               aes(y = dat_quantile$mild_mod_k1[3],
                   yend = dat_quantile$mild_mod_k1[3],
                   x = 1.75, xend = 1.85), linewidth = 0.4, inherit.aes = F) +
  geom_segment(data=data.frame(pet_param="avg_k1"), 
               aes(y = dat_quantile$severe_k1[3],
                   yend = dat_quantile$severe_k1[3],
                   x = 2.15, xend = 2.25), linewidth = 0.4, inherit.aes = F) +
  geom_segment(data=data.frame(pet_param="avg_k1"), 
               aes(y = dat_quantile$mild_mod_k1[1],
                   yend = dat_quantile$mild_mod_k1[5],
                   x = 1.8, xend = 1.8),
               linetype = "dashed", linewidth = 0.2, inherit.aes = F) +
  geom_segment(data=data.frame(pet_param="avg_k1"), 
               aes(y = dat_quantile$severe_k1[1],
                   yend = dat_quantile$severe_k1[5],
                   x = 2.2, xend = 2.2),
               linetype = "dashed", linewidth = 0.2, inherit.aes = F) +    
  new_scale_fill() +
  geom_point(aes(fill = group_w_class),
             color = "black",
             alpha = 0.7, 
             shape = 21,
             position = position_jitterdodge()) +
  facet_grid(pet_param ~ .,
             switch = "y",
             labeller = labeller(pet_param = c("avg_f" = "Avg F", "avg_k1" = "Avg K1"))) +
  scale_fill_manual(values = c("#c2dfe3", "#fcb1a6", "#fb6376")) +
  scale_color_manual(values = c("#c2dfe3", "#fff9ec",
                                "#fcb1a6", "#fb6376")) +
  labs(x = NULL,
       y = NULL,
       fill = NULL) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        legend.position = "none",
        panel.background = element_rect(fill = "#f2e9e4",
                                        color = "grey"),
        strip.background = element_rect(fill = "#f2e9e4",
                                        color = "grey"),
        strip.text = element_text(size = 10)) +
  geom_signif(annotations = "***", y_position = c(2.5), xmin = c(1), xmax = c(2))

boxplot2 <-
  ggplot(subset(boxplot_adpkd_combined_dat,
                (pet_param == "avg_k2_w_cyst"|pet_param == "avg_k2_wo_cyst") & 
                  group_w_class != "ADPKD"),
         aes(x = group,
             y = value,
             fill = group)) + 
  scale_fill_manual(values = c("#c2dfe3", "#fff9ec", 
                               "#fcb1a6","#fb6376")) +
  labs(x = NULL,
       y = NULL,
       fill = NULL) +
  geom_boxplot(outlier.shape = NA) +
  geom_rect(data=data.frame(pet_param="avg_k2_w_cyst"), 
            aes(ymin = dat_quantile$mild_mod_k2cyst[2], 
                ymax = dat_quantile$mild_mod_k2cyst[4], 
                xmin = 1.75, 
                xmax = 1.85),
            fill = "#fcb1a6", alpha = 0.9, inherit.aes = F) +
  geom_rect(data=data.frame(pet_param="avg_k2_w_cyst"), 
            aes(ymin = dat_quantile$severe_k2cyst[2], 
                ymax = dat_quantile$severe_k2cyst[4], 
                xmin = 2.15, 
                xmax = 2.25),
            fill = "#fb6376", alpha = 0.9, inherit.aes = F) +
  geom_segment(data=data.frame(pet_param="avg_k2_w_cyst"), 
               aes(y = dat_quantile$mild_mod_k2cyst[3],
                   yend = dat_quantile$mild_mod_k2cyst[3],
                   x = 1.75, xend = 1.85), linewidth = 0.4, inherit.aes = F) +
  geom_segment(data=data.frame(pet_param="avg_k2_w_cyst"), 
               aes(y = dat_quantile$severe_k2cyst[3],
                   yend = dat_quantile$severe_k2cyst[3],
                   x = 2.15, xend = 2.25), linewidth = 0.4, inherit.aes = F) +
  geom_segment(data=data.frame(pet_param="avg_k2_w_cyst"), 
               aes(y = dat_quantile$mild_mod_k2cyst[1],
                   yend = dat_quantile$mild_mod_k2cyst[5],
                   x = 1.8, xend = 1.8),
               linetype = "dashed", linewidth = 0.2, inherit.aes = F) +
  geom_segment(data=data.frame(pet_param="avg_k2_w_cyst"), 
               aes(y = dat_quantile$severe_k2cyst[1],
                   yend = dat_quantile$severe_k2cyst[5],
                   x = 2.2, xend = 2.2),
               linetype = "dashed", linewidth = 0.2, inherit.aes = F) +
  geom_rect(data=data.frame(pet_param="avg_k2_wo_cyst"), 
            aes(ymin = dat_quantile$mild_mod_k2_ncyst[2], 
                ymax = dat_quantile$mild_mod_k2_ncyst[4], 
                xmin = 1.75, 
                xmax = 1.85),
            fill = "#fcb1a6", alpha = 0.9, inherit.aes = F) +
  geom_rect(data=data.frame(pet_param="avg_k2_wo_cyst"), 
            aes(ymin = dat_quantile$severe_k2_ncyst[2], 
                ymax = dat_quantile$severe_k2_ncyst[4], 
                xmin = 2.15, 
                xmax = 2.25),
            fill = "#fb6376", alpha = 0.9, inherit.aes = F) +
  geom_segment(data=data.frame(pet_param="avg_k2_wo_cyst"), 
               aes(y = dat_quantile$mild_mod_k2_ncyst[3],
                   yend = dat_quantile$mild_mod_k2_ncyst[3],
                   x = 1.75, xend = 1.85), linewidth = 0.4, inherit.aes = F) +
  geom_segment(data=data.frame(pet_param="avg_k2_wo_cyst"), 
               aes(y = dat_quantile$severe_k2_ncyst[3],
                   yend = dat_quantile$severe_k2_ncyst[3],
                   x = 2.15, xend = 2.25), linewidth = 0.4, inherit.aes = F) +
  geom_segment(data=data.frame(pet_param="avg_k2_wo_cyst"), 
               aes(y = dat_quantile$mild_mod_k2_ncyst[1],
                   yend = dat_quantile$mild_mod_k2_ncyst[5],
                   x = 1.8, xend = 1.8),
               linetype = "dashed", linewidth = 0.2, inherit.aes = F) +
  geom_segment(data=data.frame(pet_param="avg_k2_wo_cyst"), 
               aes(y = dat_quantile$severe_k2_ncyst[1],
                   yend = dat_quantile$severe_k2_ncyst[5],
                   x = 2.2, xend = 2.2),
               linetype = "dashed", linewidth = 0.2, inherit.aes = F) +   
  new_scale_fill() +
  geom_point(aes(fill = group_w_class),
             color = "black",
             alpha = 0.7, 
             shape = 21,
             position = position_jitterdodge()) +
  facet_grid(pet_param ~ .,
             switch = "y",
             labeller = labeller(pet_param = c("avg_k2_w_cyst" = "Avg k2 w/ cyst", "avg_k2_wo_cyst" = "Avg k2 w/o cyst"))) +
  scale_fill_manual(values = c("#c2dfe3", "#fcb1a6", "#fb6376")) +
  scale_color_manual(values = c("#c2dfe3", "#fff9ec",
                                "#fcb1a6", "#fb6376")) +
  labs(x = NULL,
       y = NULL,
       fill = NULL) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        legend.position = "none",
        panel.background = element_rect(fill = "#f2e9e4",
                                        color = "grey"),
        strip.background = element_rect(fill = "#f2e9e4",
                                        color = "grey"),
        strip.text = element_text(size = 10)) +
  geom_signif(annotations = "***", y_position = c(0.205), xmin = c(1), xmax = c(2)) +
  ylim(c(0.12,0.215))

layout <- c(
  area(1, 3, 2, 6),
  area(3, 1, 6, 4),
  area(3, 5, 6, 9)
)
plot(layout)
bmi_grp + boxplot1+boxplot2 +
  plot_layout(design = layout, guides = "collect") & theme(legend.position = "top")

ggsave("/Users/choiyej/GitHub/YC_CHCO/PENGUIN/boxplots.png", width = 7, height = 7)