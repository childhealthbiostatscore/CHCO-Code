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


dat_results <- 


boxplot_function <- function(data, combined_data, variable, label){

  data <- data %>% filter()
  
  ggplot(data, aes(x = group, y = value, fill = group),
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
  
}







boxplot1 <-
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