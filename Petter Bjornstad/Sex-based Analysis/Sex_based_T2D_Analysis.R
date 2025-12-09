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




load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')

so_subset <- so_kpmp_sc
remove(so_kpmp_sc)


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


test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))

dir.results <- 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/'

so_subset <- subset(so_subset, subset = record_id != 'CRC-55')
so_subset <- subset(so_subset, subset = group == 'Type_2_Diabetes')
test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

## demographics tables 

dat_small <- dat %>% filter(record_id %in% test$record_id)

library(gt)
library(gtsummary)

desc_table1_fixed <- dat_small %>%
  select(age, sex, race_ethnicity, bmi, hba1c, study) %>%
  tbl_summary(
    by = sex,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      hba1c ~ "continuous",
      race_ethnicity ~ "categorical",
      study ~ "categorical"
      
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      study ~ "Study"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "t.test"
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Group**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

desc_table1_fixed %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave(paste0(dir.results, "T2D_demographics.png"), 
         vwidth = 1200, vheight = 800)

#### Celltype Analyses

## Barplots and UMAPs 

dir.results <- '/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/'

png(paste0(dir.results, 'T2D_UMAP.png'), width = 1200, height = 1200)
DimPlot(so_subset) + 
  xlab("UMAP 1") + 
  ylab("UMAP 2") +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  )
dev.off()

# Extract cell type and group information
cell_data <- data.frame(
  cell_type = so_subset$celltype2,
  group = so_subset$sex
)
cell_data <- cell_data %>% filter(cell_type %in% c('EC', 'PT', 'TAL'))

# Calculate proportions
prop_data <- cell_data %>%
  group_by(group, cell_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(group) %>%
  mutate(
    total = sum(count),
    proportion = count / total * 100
  )

# Statistical testing
stat_results <- data.frame()

for(ct in unique(cell_data$cell_type)) {
  contingency <- table(
    cell_data$group,
    cell_data$cell_type == ct
  )
  
  test <- chisq.test(contingency)
  
  stat_results <- rbind(stat_results, data.frame(
    cell_type = ct,
    p_value = test$p.value
  ))
}

stat_results$p_adj <- p.adjust(stat_results$p_value, method = "BH")
stat_results$sig_label <- case_when(
  stat_results$p_adj < 0.001 ~ "***",
  stat_results$p_adj < 0.01 ~ "**",
  stat_results$p_adj < 0.05 ~ "*",
  TRUE ~ ""
)

# Calculate cumulative proportions for positioning asterisks
prop_data_cumsum <- prop_data %>%
  arrange(group, desc(cell_type)) %>%
  group_by(group) %>%
  mutate(
    pos = cumsum(proportion) - proportion/2
  )

# Add significance labels to proportion data
prop_data_cumsum <- prop_data_cumsum %>%
  left_join(stat_results %>% select(cell_type, sig_label), by = "cell_type")

# Create plot with asterisks on each bar segment
p <- ggplot(prop_data, aes(x = group, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.3) +
  geom_text(data = prop_data_cumsum %>% filter(sig_label != ""),
            aes(x = group, y = pos, label = sig_label),
            size = 8, color = "white", fontface = "bold") +
  labs(
    x = "Sex in T2D",
    y = "Proportion (%)",
    fill = "Cell Type",
    title = paste0("Significant cell types (p<0.05): ", 
                   paste(stat_results$cell_type[stat_results$p_adj < 0.05], 
                         collapse = ", "))
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 11)
  )

ggsave(paste0(dir.results, 'cell_type_proportions_with_asterisks.png'), 
       plot = p, width = 8, height = 6, dpi = 300)

# Alternative: Add asterisks above the bars
prop_data_top <- prop_data %>%
  arrange(group, desc(cell_type)) %>%
  group_by(group) %>%
  mutate(y_top = cumsum(proportion)) %>%
  ungroup() %>%
  left_join(stat_results %>% select(cell_type, sig_label), by = "cell_type") %>%
  filter(sig_label != "")

p2 <- ggplot(prop_data, aes(x = group, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.3) +
  geom_text(data = prop_data_top,
            aes(x = group, y = y_top + 2, label = sig_label),
            size = 6, color = "black", fontface = "bold") +
  labs(
    x = "Sex in T2D",
    y = "Proportion (%)",
    fill = "Cell Type"
  ) +
  ylim(0, 105) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

ggsave(paste0(dir.results, 'cell_type_proportions_asterisks_above.png'), 
       plot = p2, width = 8, height = 6, dpi = 300)

print(stat_results)

for(celltype in c('TAL', 'PT', 'EC')){
  
  cell_data <- data.frame(
    overall_cell = so_subset$celltype2,
    cell_type = so_subset$KPMP_celltype,
    group = so_subset$sex
  )
  cell_data <- cell_data %>% filter(overall_cell == celltype) %>%
    dplyr::select(-overall_cell)
  
  # Calculate proportions
  prop_data <- cell_data %>%
    group_by(group, cell_type) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(group) %>%
    mutate(
      total = sum(count),
      proportion = count / total * 100
    )
  
  # Statistical testing
  stat_results <- data.frame()
  
  for(ct in unique(cell_data$cell_type)) {
    contingency <- table(
      cell_data$group,
      cell_data$cell_type == ct
    )
    
    test <- chisq.test(contingency)
    
    stat_results <- rbind(stat_results, data.frame(
      cell_type = ct,
      p_value = test$p.value
    ))
  }
  
  stat_results$p_adj <- p.adjust(stat_results$p_value, method = "BH")
  stat_results$sig_label <- case_when(
    stat_results$p_adj < 0.001 ~ "***",
    stat_results$p_adj < 0.01 ~ "**",
    stat_results$p_adj < 0.05 ~ "*",
    TRUE ~ ""
  )
  
  # Calculate cumulative proportions for positioning asterisks
  prop_data_cumsum <- prop_data %>%
    arrange(group, desc(cell_type)) %>%
    group_by(group) %>%
    mutate(
      pos = cumsum(proportion) - proportion/2
    )
  
  # Add significance labels to proportion data
  prop_data_cumsum <- prop_data_cumsum %>%
    left_join(stat_results %>% select(cell_type, sig_label), by = "cell_type")
  
  # Create plot with asterisks on each bar segment
  p <- ggplot(prop_data, aes(x = group, y = proportion, fill = cell_type)) +
    geom_bar(stat = "identity", position = "stack", color = "black", size = 0.3) +
    geom_text(data = prop_data_cumsum %>% filter(sig_label != ""),
              aes(x = group, y = pos, label = sig_label),
              size = 8, color = "white", fontface = "bold") +
    labs(
      x = "Sex in T2D",
      y = "Proportion (%)",
      fill = "Cell Type",
      title = paste0("Significant cell types (p<0.05): ", 
                     paste(stat_results$cell_type[stat_results$p_adj < 0.05], 
                           collapse = ", "))
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      plot.title = element_text(size = 11)
    )
  
  ggsave(paste0(dir.results, celltype, '_proportions_with_asterisks.png'), 
         plot = p, width = 8, height = 6, dpi = 300)
  
  # Alternative: Add asterisks above the bars
  prop_data_top <- prop_data %>%
    arrange(group, desc(cell_type)) %>%
    group_by(group) %>%
    mutate(y_top = cumsum(proportion)) %>%
    ungroup() %>%
    left_join(stat_results %>% select(cell_type, sig_label), by = "cell_type") %>%
    filter(sig_label != "")
  
  p2 <- ggplot(prop_data, aes(x = group, y = proportion, fill = cell_type)) +
    geom_bar(stat = "identity", position = "stack", color = "black", size = 0.3) +
    geom_text(data = prop_data_top,
              aes(x = group, y = y_top + 2, label = sig_label),
              size = 6, color = "black", fontface = "bold") +
    labs(
      x = "Sex in T2D",
      y = "Proportion (%)",
      fill = "Cell Type"
    ) +
    ylim(0, 105) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12)
    )
  
  ggsave(paste0(dir.results, celltype, '_proportions_asterisks_above.png'), 
         plot = p2, width = 8, height = 6, dpi = 300)
  
  print(stat_results)
  print(celltype)
}

################## aPT analysis

meta.data <- so_subset@meta.data

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

png(paste0(dir.results, "aPT_Percentage_byParticipant.pdf"))
ggplot(cell_distribution, aes(x=group_labels, color=group_labels, y= aPT_percentage))+
  geom_point(position=position_jitter(width = 0.1, height = 0))+geom_boxplot()+theme_classic()+labs(x='Condition Group', y='Percent aPT of PT Cells', title = 'aPT Percentage in Each Participant by Group')
dev.off()

png(paste0(dir.results, "aPT_CountbyParticipant.png"))
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

png(paste0(dir.results, "dTAL_Percentage_byParticipant.png"))
ggplot(cell_distribution, aes(x=group_labels, color=group_labels, y= dTAL_percentage))+
  geom_point(position=position_jitter(width = 0.1, height = 0))+geom_boxplot()+theme_classic()+labs(x='Condition Group', y='Percent dTAL in TAL Cells', title = 'dTAL Percentage in Each Participant by Group')
dev.off()

png(paste0(dir.results, "dTAL_CountbyParticipant.png"))
ggplot(cell_distribution, aes(x=group_labels, color=group_labels, y= dTAL_count))+
  geom_point(position=position_jitter(width = 0.1, height = 0))+geom_boxplot()+theme_classic()+labs(x='Condition Group', y='Number of dTAL Cells', title = 'dTAL Cells in Each Participant by Group')
dev.off()

TAL_data <- cell_distribution

#Combined for analyses

names(aPT_data) <- c('record_id', 'group_labels', 'PT_totalcells', 'aPT_count', 'aPT_percentage')
names(TAL_data) <- c('record_id', 'group_labels', 'TAL_totalcells', 'dTAL_count', 'dTAL_percentage')

strange_cells <- aPT_data %>% left_join(TAL_data)
sex_df <- meta.data %>% dplyr::select(record_id, sex)

strange_cells <- strange_cells %>% left_join(sex_df, by='record_id')

ggplot(strange_cells, aes(x=aPT_percentage, y=dTAL_percentage, color = sex))+geom_point()+
  geom_smooth(method = 'lm')+theme_classic()+theme(x = 'aPT Percentage of PT Cells', 
                                                   y = 'dTAL Percentage of TAL Cells')

lm(dTAL_percentage ~ aPT_percentage*sex, data = strange_cells) %>% summary()

library(ggplot2)
library(broom)

# Model
model <- lm(dTAL_percentage ~ aPT_percentage * sex, data = strange_cells)
model_summary <- tidy(model)

# Extract key statistics
interaction_p <- model_summary$p.value[model_summary$term == "aPT_percentage:sexMale"]

# Significance star
sig_star <- case_when(
  interaction_p < 0.001 ~ "***",
  interaction_p < 0.01 ~ "**",
  interaction_p < 0.05 ~ "*",
  TRUE ~ "ns"
)

# Plot
p <- ggplot(strange_cells, aes(x = aPT_percentage, y = dTAL_percentage, 
                               color = sex, fill = sex)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = 'lm', se = TRUE, alpha = 0.15, linewidth = 1.3) +
  scale_color_manual(values = c("Male" = "#0073C2FF", "Female" = "#EFC000FF")) +
  scale_fill_manual(values = c("Male" = "#0073C2FF", "Female" = "#EFC000FF")) +
  labs(
    x = 'aPT Percentage of PT Cells (%)', 
    y = 'dTAL Percentage of TAL Cells (%)',
    color = 'Sex',
    fill = 'Sex',
    title = "Sex-specific relationship between aPT and dTAL cell populations"
  ) +
  annotate("text", x = Inf, y = Inf, 
           label = sprintf("Interaction: p = %.4f %s", interaction_p, sig_star),
           hjust = 1.1, vjust = 2, size = 4, fontface = "bold") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.position = c(0.85, 0.15),
    legend.background = element_rect(fill = "white", color = "black")
  )

print(p)
print(summary(model))

ggsave(paste0(dir.results, "aPT_dTAL_interaction_plot.png"), plot = p, 
       width = 8, height = 6, dpi = 300)

pathway_data <- data.table::fread("/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/omics/metabolomics_pathways/pathway_results.csv")

# Filter for significant pathways
significant_pathways <- pathway_data[pathway_data$`Raw p` < 0.05, ]

# Create shortened labels
shorten_label <- function(label, max_length = 30) {
  ifelse(nchar(label) > max_length, 
         paste0(substr(label, 1, max_length), "..."), 
         label)
}

significant_pathways$label_short <- shorten_label(significant_pathways$V1)

png('/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/omics/metabolomics_pathways/metabolomics_pathways.png', 
    width = 1200, height = 1200)

ggplot(pathway_data, aes(x = Impact, y = `-log10(p)`)) +
  geom_point(aes(size = Impact, fill = `-log10(p)`), 
             shape = 21, 
             color = "black", 
             stroke = 0.5) +
  geom_text_repel(data = significant_pathways,
                  aes(label = label_short),
                  size = 6,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  max.overlaps = Inf,
                  segment.color = "grey50") +
  scale_size_continuous(range = c(2, 15)) +
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(x = "Pathway Impact", 
       y = "-log10(p)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

dev.off()

#function
T2D_NEBULA_Analysis <- function(so_subset, dir.results, celltype, genes){
  if(celltype == 'All'){
    so_celltype <- so_subset 
  }else if(celltype %in% c('TAL', 'EC', 'POD', 'PT')){
    so_celltype <- subset(so_subset,celltype2==celltype)
    DefaultAssay(so_celltype) <- "RNA" 
  }else if(celltype == 'IC'){
    so_celltype <- subset(so_subset, KPMP_celltype %in% c(
      "cDC",
      "cycT",
      "CD4+ T",
      "CD8+ T",
      "NK",
      "B",
      "MON",
      "MAC",
      "MC"))
  }else if(celltype == 'DCTall'){
    so_celltype <- subset(so_subset, DCT_celltype=='DCT')
  }else if(celltype == 'intercalated'){
    so_celltype <- subset(so_subset, KPMP_celltype %in% c('IC-B', 'IC-A', 'tPC-IC', 'aIC'))
  }
  
  full_analysis <- FindVariableFeatures(so_celltype, selection.method = "vst", nfeatures = 2000)
  hvgs <- VariableFeatures(full_analysis)
  so_celltype <- subset(so_celltype, features = hvgs)
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  print(paste0(celltype2, ' is running.'))
  
  counts_path <- round(GetAssayData(so_celltype, layer = "counts"))
  count_gene <- counts_path
  meta_gene <- subset(so_celltype)@meta.data
  
  complete_idx <- complete.cases(meta_gene$sex)
  cat("Cells with complete data:", sum(complete_idx), "\n")
  
  meta_gene <- meta_gene[complete_idx, ]
  count_gene <- count_gene[, complete_idx]
  
  num_cells <- nrow(meta_gene)
  num_part <- unique(meta_gene$record_id) %>% length()
  
  tmp_df <- meta_gene %>% dplyr::select(record_id, group, sex) %>% filter(!duplicated(record_id))
  
  num_male <- tmp_df %>% filter(sex == 'Male') %>% nrow()
  num_female <- tmp_df %>% filter(sex == 'Female') %>% nrow()
  
  pred_gene <- model.matrix(~sex, data = meta_gene)
  
  library <- meta_gene$pooled_offset
  data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)
  
  if (is.null(data_g_gene)) {
    data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
  }
  
  result <- nebula(count = data_g_gene$count, id = data_g_gene$id, 
                   pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)
  
  full_results <- as.data.frame(result)
  
  full_results$num_cells <- num_cells
  full_results$num_male <- num_male
  full_results$num_female <- num_female
  
  write.table(full_results,paste0(dir.results,"NEBULA_", 
                                  celltype2, "_cells__T2D_pooledoffset.csv"),
              row.names=F, quote=F, sep=',')
  
  print(paste0(celltype2, ' is done.'))
}

celltypes_vec <- c('All', 
                   'PT', 
                   'TAL', 
                   'EC', 
                   'POD',
                   'DCTall', 
                   'IC')

for(i  in 1:length(celltypes_vec)){
  T2D_NEBULA_Analysis(so_subset = so_subset, dir.results = dir.results, celltype = celltypes_vec[i], genes = gene_list)
  print(celltypes_vec[i])
}

### Full Analysis (all Genes)

T2D_NEBULA_Analysis_full <- function(so_subset, dir.results, celltype, genes){
  if(celltype == 'All'){
    so_celltype <- so_subset 
  }else if(celltype %in% c('TAL', 'EC', 'POD', 'PT')){
    so_celltype <- subset(so_subset,celltype2==celltype)
    DefaultAssay(so_celltype) <- "RNA" 
  }else if(celltype == 'IC'){
    so_celltype <- subset(so_subset, KPMP_celltype %in% c(
      "cDC",
      "cycT",
      "CD4+ T",
      "CD8+ T",
      "NK",
      "B",
      "MON",
      "MAC",
      "MC"))
  }else if(celltype == 'DCTall'){
    so_celltype <- subset(so_subset, DCT_celltype=='DCT')
  }else if(celltype == 'intercalated'){
    so_celltype <- subset(so_subset, KPMP_celltype %in% c('IC-B', 'IC-A', 'tPC-IC', 'aIC'))
  }
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  print(paste0(celltype2, ' is running.'))
  
  counts_path <- round(GetAssayData(so_celltype, layer = "counts"))
  count_gene <- counts_path
  meta_gene <- subset(so_celltype)@meta.data
  
  complete_idx <- complete.cases(meta_gene$sex)
  cat("Cells with complete data:", sum(complete_idx), "\n")
  
  meta_gene <- meta_gene[complete_idx, ]
  count_gene <- count_gene[, complete_idx]
  
  num_cells <- nrow(meta_gene)
  num_part <- unique(meta_gene$record_id) %>% length()
  
  tmp_df <- meta_gene %>% dplyr::select(record_id, group, sex) %>% filter(!duplicated(record_id))
  
  num_male <- tmp_df %>% filter(sex == 'Male') %>% nrow()
  num_female <- tmp_df %>% filter(sex == 'Female') %>% nrow()
  
  pred_gene <- model.matrix(~sex, data = meta_gene)
  
  library <- meta_gene$pooled_offset
  data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)
  
  if (is.null(data_g_gene)) {
    data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
  }
  
  result <- nebula(count = data_g_gene$count, id = data_g_gene$id, 
                   pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)
  
  full_results <- as.data.frame(result)
  
  full_results$num_cells <- num_cells
  full_results$num_male <- num_male
  full_results$num_female <- num_female
  
  write.table(full_results,paste0(dir.results,"Full_NEBULA_", 
                                  celltype2, "_cells__T2D_pooledoffset.csv"),
              row.names=F, quote=F, sep=',')
  
  print(paste0(celltype2, ' is done.'))
}

celltypes_vec <- c('All', 
                   'PT', 
                   'TAL', 
                   'EC', 
                   'POD',
                   'DCTall', 
                   'IC')

for(i  in 1:length(celltypes_vec)){
  T2D_NEBULA_Analysis_full(so_subset = so_subset, dir.results = dir.results, celltype = celltypes_vec[i], genes = gene_list)
  print(celltypes_vec[i])
}

########### Volcano plots 

library(ggplot2)
library(ggbreak)
library(dplyr)

variable_names <- c('All', 'PT', 'TAL', 'EC', 'IC', 'POD', 'DCTall', 'intercalated')

for(i in c(1:length(variable_names))){
  
  sig_markers <- data.table::fread(paste0('/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/Full_NEBULA_', 
                                          variable_names[i], '_cells__T2D_pooledoffset.csv'))
  
  sig_markers <- sig_markers %>% dplyr::select(Gene = summary.gene,
                                               LogFC = summary.logFC_sexMale, 
                                               Pvalue = summary.p_sexMale)
  
  tmp_df <- sig_markers
  tmp_df$diffexp <- 'No'
  tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC > 0] <- 'Up'
  tmp_df$diffexp[tmp_df$Pvalue < 0.05 & tmp_df$LogFC < 0] <- 'Down'
  
  tmp_df <- tmp_df %>% arrange(Pvalue)
  tmp_df$label <- NA
  tmp_df$label[1:10] <- tmp_df$Gene[1:10]
  
  tmp_df <- tmp_df %>% filter(abs(LogFC) < 10)
  
  logfc_range <- range(tmp_df$LogFC, na.rm = TRUE)
  logfc_q95 <- quantile(abs(tmp_df$LogFC), 0.95, na.rm = TRUE)
  logfc_max <- max(abs(tmp_df$LogFC), na.rm = TRUE)
  
  needs_break <- logfc_max > (logfc_q95 * 2)
  
  if(needs_break){
    sorted_logfc <- sort(abs(tmp_df$LogFC))
    gaps <- diff(sorted_logfc)
    gap_threshold <- quantile(abs(tmp_df$LogFC), 0.90, na.rm = TRUE)
    outer_data_idx <- which(sorted_logfc > gap_threshold)
    
    if(length(outer_data_idx) > 1){
      outer_gaps <- gaps[outer_data_idx[-length(outer_data_idx)]]
      max_gap_idx <- which.max(outer_gaps)
      
      break_start <- sorted_logfc[outer_data_idx[max_gap_idx]] + 0.1
      break_end <- sorted_logfc[outer_data_idx[max_gap_idx + 1]] - 0.1
      
      if(break_end - break_start > 0.5){
        use_break <- TRUE
      } else {
        use_break <- FALSE
      }
    } else {
      use_break <- FALSE
    }
  } else {
    use_break <- FALSE
  }
  
  if(length(unique(tmp_df$diffexp)) > 1){
    tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('orange', 'grey', 'purple'),
                         labels = c('Downregulated', 'Not significant', 'Upregulated'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
      theme_classic()+
      labs(x='LogFC', y='-log10 pvalue', col ='Differential Expression', 
           title = paste0('Sex Differences in ',  variable_names[i],' Cells'))
    
    if(use_break){
      tmp_graph <- tmp_graph + scale_x_break(breaks = c(break_start, break_end), scales = 0.3)
    }
    
  } else {
    tmp_graph <- ggplot(tmp_df, aes(x= LogFC, y=-log10(Pvalue), col = diffexp, label=label))+
      geom_point()+
      geom_text(size=2, vjust = 2, color='black')+
      scale_color_manual(values = c('grey'),
                         labels = c('Not significant'))+
      geom_hline(yintercept = -log10(0.05), col='blue', linetype='dashed')+
      geom_vline(xintercept = c(0), col='black', linetype ='dashed')+
      theme_classic()+
      labs(x='LogFC', y='-log10 P-value', col ='Differential Expression', 
           title = paste0('Sex Differences in ',  variable_names[i],' Cells'))
    
    if(use_break){
      tmp_graph <- tmp_graph + scale_x_break(breaks = c(break_start, break_end), scales = 0.3)
    }
  }
  
  pdf(paste0('/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/VolcanoPlots_', 
             variable_names[i], '_Cells_T2DOnly.pdf'))
  print(tmp_graph)
  dev.off()
  
  png(paste0('/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/VolcanoPlots_', 
             variable_names[i], '_Cells_T2DOnly.png'), width = 3000, height = 3000, res = 300)
  print(tmp_graph)
  dev.off()
  
  print(paste0('Plot done for ', variable_names[i], 
               if(use_break) ' (with axis break)' else ''))
}

##### GSEA

library(fgsea)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(patchwork)

dir.results <- 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/'
celltypes_vec <- c('All', 
                   'PT', 
                   'TAL', 
                   'EC', 
                   'POD',
                   'DCTall', 
                   'IC')

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

go_bp_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
go_bp_list <- split(go_bp_sets$gene_symbol, go_bp_sets$gs_name)

go_cc_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")
go_cc_list <- split(go_cc_sets$gene_symbol, go_cc_sets$gs_name)

go_mf_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")
go_mf_list <- split(go_mf_sets$gene_symbol, go_mf_sets$gs_name)

geneset_types <- list(
  Hallmark = hallmark_list,
  GO_BP = go_bp_list,
  GO_CC = go_cc_list,
  GO_MF = go_mf_list
)

all_gsea_results <- list()

for(celltype in celltypes_vec) {
  
  cat("\n=== Processing", celltype, "===\n")
  
  de_results <- read.csv(paste0(dir.results, 'Full_NEBULA_', 
                                celltype, '_cells__T2D_pooledoffset.csv'))
  
  ranked_genes <- de_results %>%
    dplyr::select(logFC = summary.logFC_sexMale, 
                  gene_id = summary.gene) %>% 
    arrange(desc(logFC)) %>%
    pull(logFC, name = gene_id)
  
  for(geneset_name in names(geneset_types)) {
    
    cat("Running GSEA for", geneset_name, "...\n")
    
    gsea_results <- fgsea(
      pathways = geneset_types[[geneset_name]],
      stats = ranked_genes,
      minSize = 15,
      maxSize = 500,
      nperm = 10000
    )
    
    gsea_results$celltype <- celltype
    gsea_results$geneset_type <- geneset_name
    
    result_name <- paste0(celltype, "_", geneset_name)
    all_gsea_results[[result_name]] <- gsea_results
    
    top_pathways <- gsea_results %>%
      filter(pval < 0.05) %>%
      arrange(pval) %>%
      slice_head(n = 20) %>%
      mutate(pathway_clean = gsub("HALLMARK_", "", pathway),
             pathway_clean = gsub("GOBP_", "", pathway_clean),
             pathway_clean = gsub("GOCC_", "", pathway_clean),
             pathway_clean = gsub("GOMF_", "", pathway_clean),
             pathway_clean = gsub("_", " ", pathway_clean),
             pathway_clean = tools::toTitleCase(tolower(pathway_clean)))
    
    if(nrow(top_pathways) > 0) {
      
      p <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway_clean, NES))) +
        geom_point(aes(size = -log10(pval), color = NES)) +
        scale_color_gradient2(
          low = "blue", 
          mid = "white", 
          high = "red", 
          midpoint = 0, 
          name = "NES"
        ) +
        scale_size_continuous(
          name = "-log10(pval)",
          range = c(2, 10)
        ) +
        theme_bw() +
        theme(
          axis.text.y = element_text(size = 9, color = "black"),
          axis.text.x = element_text(size = 9, color = "black"),
          panel.grid.major.y = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
        ) +
        labs(
          x = "Normalized Enrichment Score (NES)",
          y = "",
          title = paste0(celltype, " - ", geneset_name, " (Top 20)")
        ) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5)
      
      ggsave(paste0(dir.results, 'GSEA/', geneset_name, '_gsea_', celltype, ".pdf"), 
             plot = p, width = 10, height = 8)
      ggsave(paste0(dir.results, 'GSEA/', geneset_name, '_gsea_', celltype, ".png"), 
             plot = p, width = 10, height = 8, dpi = 300)
      
      cat("Saved plot for", celltype, "-", geneset_name, "\n")
    } else {
      cat("No significant pathways found for", celltype, "-", geneset_name, "\n")
    }
  }
}

combined_gsea <- bind_rows(all_gsea_results)

combined_gsea_save <- combined_gsea %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";")))

write.csv(combined_gsea_save, 
          paste0(dir.results, "GSEA/all_celltypes_all_genesets_gsea_results.csv"), 
          row.names = FALSE)

for(geneset_name in names(geneset_types)) {
  subset_results <- combined_gsea %>% 
    filter(geneset_type == geneset_name) %>%
    mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";")))
  
  write.csv(subset_results, 
            paste0(dir.results, "GSEA/gsea_results_", geneset_name, ".csv"), 
            row.names = FALSE)
}

cat("\n=== All analyses complete! ===\n")
cat("Total number of analyses:", length(all_gsea_results), "\n")

summary_table <- combined_gsea %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";"))) %>% 
  filter(pval < 0.05) %>%
  group_by(celltype, geneset_type) %>%
  summarise(n_significant = n(), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = geneset_type, values_from = n_significant, values_fill = 0)

print(summary_table)
write.csv(summary_table, 
          paste0(dir.results, "GSEA/summary_significant_pathways.csv"), 
          row.names = FALSE)

# Convert PDFs to PNGs
if (!require("pdftools")) {
  install.packages("pdftools")
}

library(pdftools)

convert_pdf_to_png <- function(pdf_path, output_dir = NULL, dpi = 300) {
  if (is.null(output_dir)) {
    output_dir <- dirname(pdf_path)
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  base_name <- tools::file_path_sans_ext(basename(pdf_path))
  
  png_files <- pdf_convert(
    pdf = pdf_path,
    format = "png",
    dpi = dpi,
    filenames = file.path(output_dir, paste0(base_name, "_page_%d.png"))
  )
  
  return(png_files)
}

convert_all_pdfs <- function(folder_path, output_dir = NULL, dpi = 300) {
  pdf_files <- list.files(folder_path, pattern = "\\.pdf$", 
                          full.names = TRUE, ignore.case = TRUE)
  
  if (length(pdf_files) == 0) {
    message("No PDF files found in the specified folder.")
    return(NULL)
  }
  
  message(paste("Found", length(pdf_files), "PDF file(s)"))
  
  all_png_files <- list()
  for (pdf in pdf_files) {
    message(paste("Converting:", basename(pdf)))
    png_files <- convert_pdf_to_png(pdf, output_dir, dpi)
    all_png_files[[basename(pdf)]] <- png_files
  }
  
  message("Conversion complete!")
  return(all_png_files)
}

folder_path <- "/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/GSEA/"
converted_files <- convert_all_pdfs(folder_path, output_dir = "/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/GSEA/")

#### Proteomics and Metabolomics 

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

dat <- dat %>% filter(group == 'Type 2 Diabetes') %>% 
  filter(visit == 'baseline')

data_dictionary <- readxl::read_xlsx('/Users/netio/Downloads/data_dictionary_master.xlsx')

form_names <- unique(data_dictionary$form_name)
proteo <- form_names[str_which(form_names, pattern = 'proteom')]
metab <- form_names[str_which(form_names, pattern = 'metab')]

variables_class <- c('proteomics', 'metabolomics', 'metabolomics_blood_raw',
                     'az_urine_metabolites', 'metabolomics_aq')

data_dictionary_small <- data_dictionary %>% 
  filter(form_name %in% variables_class)

library(gtsummary)
library(gt)
library(dplyr)

combined_df <- dat %>% 
  filter(!is.na(citrate)) %>% 
  mutate(
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    hba1c = as.numeric(hba1c),
    sex = as.factor(sex),
    race_ethnicity = as.factor(race_ethnicity),
    study = as.factor(study),
    group = as.factor(group)
  )

desc_table1_fixed <- combined_df %>%
  select(age, sex, race_ethnicity, bmi, hba1c, study, group) %>%
  tbl_summary(
    by = sex,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      hba1c ~ "continuous",
      group ~ "categorical",
      race_ethnicity ~ "categorical",
      study ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      group ~ "Group", 
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      study ~ "Study"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "t.test"
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Sex**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

desc_table1_fixed %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave("C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/omics/T2D_Omics_Demographics.png", 
         vwidth = 1200, vheight = 800)

existing_cols <- intersect(data_dictionary_small$variable_name, names(dat))

dat_omics <- dat %>% 
  dplyr::select(record_id, group, sex, study, all_of(existing_cols))

cat("Found:", length(existing_cols), "out of", length(data_dictionary_small$variable_name), "\n")
cat("Missing:", length(data_dictionary_small$variable_name) - length(existing_cols), "\n")

##Analysis

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(readxl)

data_dictionary <- readxl::read_xlsx('/Users/netio/Downloads/data_dictionary_master.xlsx')

output_folder <- '/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/omics/'

numeric_cols <- names(dat_omics)[!names(dat_omics) %in% c('record_id', 'group', 'sex')]

results <- data.frame(
  variable = character(),
  variable_label = character(),
  form_name = character(),
  test_used = character(),
  p_value = numeric(),
  mean_male = numeric(),
  mean_female = numeric(),
  stringsAsFactors = FALSE
)

for(var in numeric_cols){
  
  var_label <- data_dictionary$label[data_dictionary$variable_name == var]
  var_form <- data_dictionary$form_name[data_dictionary$variable_name == var]
  
  if(length(var_label) == 0) var_label <- var
  if(length(var_form) == 0) var_form <- NA
  
  temp_data <- dat_omics %>%
    filter(!is.na(sex) & sex != "") %>%
    select(sex, all_of(var)) %>%
    rename(value = all_of(var)) %>%
    filter(!is.na(value))
  
  if(nrow(temp_data) < 6 || length(unique(temp_data$sex)) < 2) next
  
  male_data <- temp_data$value[temp_data$sex == "Male"]
  female_data <- temp_data$value[temp_data$sex == "Female"]
  
  if(length(male_data) < 3 || length(female_data) < 3) next
  
  normal_male <- if(length(male_data) >= 3 & length(male_data) <= 5000) {
    shapiro.test(male_data)$p.value > 0.05
  } else TRUE
  
  normal_female <- if(length(female_data) >= 3 & length(female_data) <= 5000) {
    shapiro.test(female_data)$p.value > 0.05
  } else TRUE
  
  if(normal_male & normal_female){
    test_result <- t.test(value ~ sex, data = temp_data, var.equal = FALSE)
    test_name <- "Welch's t-test"
  } else {
    test_result <- wilcox.test(value ~ sex, data = temp_data)
    test_name <- "Mann-Whitney U"
  }
  
  results <- rbind(results, data.frame(
    variable = var,
    variable_label = var_label,
    form_name = var_form,
    test_used = test_name,
    p_value = test_result$p.value,
    mean_male = mean(male_data, na.rm = TRUE),
    mean_female = mean(female_data, na.rm = TRUE)
  ))
}

results <- results %>%
  mutate(
    p_bonferroni = p.adjust(p_value, method = "bonferroni"),
    p_fdr = p.adjust(p_value, method = "fdr"),
    difference = mean_male - mean_female,
    abs_difference = abs(difference)
  ) %>%
  arrange(p_value)

print(head(results, 20))

write.csv(results, paste0(output_folder, 'sex_comparison_results.csv'), row.names = FALSE)

top_results <- head(results, 10)

plot_list <- list()

for(i in 1:nrow(top_results)){
  var <- top_results$variable[i]
  var_label <- top_results$variable_label[i]
  
  plot_data <- dat_omics %>%
    filter(!is.na(sex) & sex != "") %>%
    select(sex, all_of(var)) %>%
    rename(value = all_of(var)) %>%
    filter(!is.na(value))
  
  p_val <- top_results$p_value[i]
  p_fdr <- top_results$p_fdr[i]
  test_type <- top_results$test_used[i]
  
  p_label <- ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val))
  p_fdr_label <- ifelse(p_fdr < 0.001, "FDR < 0.001", sprintf("FDR = %.3f", p_fdr))
  
  plot_list[[i]] <- ggplot(plot_data, aes(x = sex, y = value, fill = sex)) +
    geom_boxplot() +
    scale_fill_manual(values = c("Male" = "#00008B", "Female" = "#8B0000")) +
    theme_classic(base_size = 12) +
    labs(
      title = var_label,
      subtitle = paste0(p_label, " | ", p_fdr_label),
      x = "Sex",
      y = "Value",
      fill = "Sex"
    ) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8),
      legend.position = "none"
    )
}

png(paste0(output_folder, 'top10_sex_differences.png'), 
    height = 15, width = 12, units = 'in', res = 300)
grid.arrange(grobs = plot_list, ncol = 2)
dev.off()

top_results_20 <- head(results, 20)
plot_list_20 <- list()

for(i in 1:nrow(top_results_20)){
  var <- top_results_20$variable[i]
  var_label <- top_results_20$variable_label[i]
  
  plot_data <- dat_omics %>%
    filter(!is.na(sex) & sex != "") %>%
    select(sex, all_of(var)) %>%
    rename(value = all_of(var)) %>%
    filter(!is.na(value))
  
  p_val <- top_results_20$p_value[i]
  p_fdr <- top_results_20$p_fdr[i]
  
  p_label <- ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val))
  p_fdr_label <- ifelse(p_fdr < 0.001, "FDR < 0.001", sprintf("FDR = %.3f", p_fdr))
  
  plot_list_20[[i]] <- ggplot(plot_data, aes(x = sex, y = value, fill = sex)) +
    geom_boxplot() +
    scale_fill_manual(values = c("Male" = "#00008B", "Female" = "#8B0000")) +
    theme_classic(base_size = 10) +
    labs(
      title = var_label,
      subtitle = paste0(p_label, " | ", p_fdr_label),
      x = "Sex",
      y = "Value",
      fill = "Sex"
    ) +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      plot.subtitle = element_text(size = 7),
      legend.position = "none"
    )
}

png(paste0(output_folder, 'top20_sex_differences.png'), 
    height = 25, width = 15, units = 'in', res = 300)
grid.arrange(grobs = plot_list_20, ncol = 3)
dev.off()

print("Analysis complete!")
print(paste("Total variables tested:", nrow(results)))
print(paste("Significant at p < 0.05:", sum(results$p_value < 0.05)))
print(paste("Significant at FDR < 0.05:", sum(results$p_fdr < 0.05)))

form_summary <- results %>%
  group_by(form_name) %>%
  summarise(
    total_variables = n(),
    sig_p05 = sum(p_value < 0.05),
    sig_fdr05 = sum(p_fdr < 0.05)
  ) %>%
  arrange(desc(sig_fdr05))
print(form_summary)




##################################### Metabolomics pathways analysis 

# Install required packages (if not already installed)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("MetaboAnalystR", "fgsea", "clusterProfiler"))
install.packages(c("FELLA", "pathview", "ggplot2", "dplyr"))

# Load libraries
library(MetaboAnalystR)
library(dplyr)
library(ggplot2)

# ===== METHOD 1: Using MetaboAnalystR (Recommended) =====

metabolites_data <- data.table::fread('/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/omics/sex_comparison_results.csv')

significant_metabolites <- metabolites_data %>%
  filter(p_value < 0.05) %>%
  pull(metabolite_name)

# Initialize MetaboAnalyst
mSet <- InitDataObjects("conc", "msetora", FALSE)

# Set up metabolite list
mSet <- Setup.MapData(mSet, significant_metabolites)

# Perform compound name mapping (KEGG IDs)
mSet <- CrossReferencing(mSet, "name")
mSet <- CreateMappingResultTable(mSet)

# Perform pathway enrichment
mSet <- SetMetabolomeFilter(mSet, F)
mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway", 2)
mSet <- CalculateHyperScore(mSet)

# Get results
pathway_results <- mSet$analSet$ora.mat
pathway_results <- as.data.frame(pathway_results)
pathway_results$pathway <- rownames(pathway_results)

# View top pathways
head(pathway_results[order(pathway_results$Raw.p),], 20)

# Visualize
ggplot(head(pathway_results[order(pathway_results$Raw.p),], 15), 
       aes(x = reorder(pathway, -log10(Raw.p)), y = -log10(Raw.p))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "Pathway", y = "-log10(p-value)", 
       title = "Top Enriched Metabolic Pathways") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))

ggsave(paste0(dir.results, "pathway_enrichment_ORA.png"), 
       width = 10, height = 8, dpi = 300)

# ===== METHOD 2: Quantitative Enrichment Analysis (MSEA) =====

metabolite_scores <- metabolites_data %>%
  select(metabolite_name, fold_change) %>%
  na.omit()

write.table(metabolite_scores, 
            file = paste0(dir.results, "metabolite_scores.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Initialize for MSEA
mSet <- InitDataObjects("conc", "msetqea", FALSE)
mSet <- Read.TextData(mSet, paste0(dir.results, "metabolite_scores.txt"), 
                      "rowu", "disc")
mSet <- SanityCheckData(mSet)
mSet <- CrossReferencing(mSet, "name")
mSet <- CreateMappingResultTable(mSet)

# Perform MSEA
mSet <- SetMetabolomeFilter(mSet, F)
mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway", 2)
mSet <- CalculateGlobalTestScore(mSet)

# Get results
msea_results <- mSet$analSet$qea.mat
msea_results <- as.data.frame(msea_results)
msea_results$pathway <- rownames(msea_results)

# ===== METHOD 3: Using FGSEA =====
library(fgsea)

metabolite_ranks <- metabolites_data %>%
  arrange(desc(fold_change)) %>%
  select(metabolite_id, fold_change) %>%
  deframe()

# Run FGSEA
fgsea_results <- fgsea(
  pathways = kegg_pathways,
  stats = metabolite_ranks,
  minSize = 5,
  maxSize = 500,
  nperm = 10000
)

fgsea_results <- fgsea_results %>%
  arrange(padj)

print(head(fgsea_results, 20))

plotEnrichment(kegg_pathways[["Glycolysis"]], metabolite_ranks) +
  labs(title = "Glycolysis Pathway Enrichment")

# ===== METHOD 4: Manual pathway enrichment using hypergeometric test =====

pathway_db <- list(
  "Glycolysis / Gluconeogenesis" = c("Glucose", "Pyruvate", "Lactate", "Glucose-6-phosphate"),
  "TCA Cycle" = c("Citrate", "Succinate", "Fumarate", "Malate"),
  "Fatty Acid Metabolism" = c("Palmitate", "Stearate", "Oleate", "Acetyl-CoA")
)

sig_metabolites <- metabolites_data %>%
  filter(p_value < 0.05) %>%
  pull(metabolite_name)

total_metabolites <- nrow(metabolites_data)

enrichment_results <- data.frame()

for(pathway_name in names(pathway_db)) {
  pathway_metabolites <- pathway_db[[pathway_name]]
  
  overlap <- length(intersect(sig_metabolites, pathway_metabolites))
  
  p_val <- phyper(
    q = overlap - 1,
    m = length(pathway_metabolites),
    n = total_metabolites - length(pathway_metabolites),
    k = length(sig_metabolites),
    lower.tail = FALSE
  )
  
  enrichment_results <- rbind(enrichment_results, data.frame(
    pathway = pathway_name,
    overlap = overlap,
    pathway_size = length(pathway_metabolites),
    p_value = p_val
  ))
}

enrichment_results$p_adj <- p.adjust(enrichment_results$p_value, method = "BH")
enrichment_results <- enrichment_results %>% arrange(p_adj)

print(enrichment_results)

# Visualize
ggplot(enrichment_results %>% filter(p_adj < 0.05), 
       aes(x = reorder(pathway, -log10(p_adj)), y = -log10(p_adj), 
           size = overlap, fill = -log10(p_adj))) +
  geom_point(shape = 21, color = "black") +
  coord_flip() +
  scale_fill_gradient(low = "lightblue", high = "darkred") +
  labs(x = "Pathway", y = "-log10(adjusted p-value)", 
       title = "Enriched Metabolic Pathways",
       size = "Number of\nMetabolites") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 11))

ggsave(paste0(dir.results, "pathway_enrichment_manual.png"), 
       width = 10, height = 8, dpi = 300)

# ===== METHOD 5: Using online MetaboAnalyst =====

write.table(significant_metabolites,
            file = paste0(dir.results, "significant_metabolites_for_metaboanalyst.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

msea_input <- metabolites_data %>%
  select(metabolite_name, fold_change) %>%
  na.omit()

write.table(msea_input,
            file = paste0(dir.results, "metabolites_with_scores.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

############## Proteomics Analyses 

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(readxl)
library(stringr)
library(httr)

# Function to get protein names from UniProt IDs
get_protein_names <- function(uniprot_ids) {
  
  unique_ids <- unique(uniprot_ids[!is.na(uniprot_ids)])
  
  if(length(unique_ids) == 0) {
    return(data.frame(Entry = character(), 
                      Gene.Names = character(), 
                      Protein.names = character()))
  }
  
  batch_size <- 100
  batches <- split(unique_ids, ceiling(seq_along(unique_ids) / batch_size))
  
  all_results <- data.frame()
  
  cat("Retrieving protein names from UniProt...\n")
  
  for (i in seq_along(batches)) {
    cat(paste0("Processing batch ", i, " of ", length(batches), "...\n"))
    
    batch <- batches[[i]]
    
    url <- "https://rest.uniprot.org/uniprotkb/search"
    query <- paste0("accession:(", paste(batch, collapse = " OR "), ")")
    
    response <- GET(
      url,
      query = list(
        query = query,
        format = "tsv",
        fields = "accession,gene_names,protein_name"
      )
    )
    
    if (status_code(response) == 200) {
      content <- content(response, "text", encoding = "UTF-8")
      if(nchar(content) > 0) {
        batch_results <- read.delim(text = content, sep = "\t", stringsAsFactors = FALSE)
        all_results <- rbind(all_results, batch_results)
      }
    } else {
      warning(paste("Failed to retrieve batch", i, ". Status code:", status_code(response)))
    }
    
    Sys.sleep(0.5)
  }
  
  cat(paste0("Retrieved information for ", nrow(all_results), " proteins\n"))
  
  return(all_results)
}

# Read harmonized data
harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/soma_harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

dat <- dat %>% filter(group == 'Type 2 Diabetes')

data_dictionary <- readxl::read_xlsx('/Users/netio/Downloads/data_dictionary_master.xlsx')

form_names <- unique(data_dictionary$form_name)
proteo <- form_names[str_which(form_names, pattern = 'proteom')]
metab <- form_names[str_which(form_names, pattern = 'metab')]

variables_class <- c('proteomics', 'az_urine_metabolites', 'metabolomics_aq')

data_dictionary_small <- data_dictionary %>% 
  filter(form_name %in% variables_class)

# Extract UniProt IDs from proteomics variable labels
cat("\n=== Extracting UniProt IDs from proteomics data ===\n")
proteomics_vars <- data_dictionary %>%
  filter(form_name == 'proteomics')

uniprot_ids <- proteomics_vars$label %>%
  str_extract("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")

protein_mapping <- get_protein_names(uniprot_ids)

data_dictionary_enhanced <- data_dictionary %>%
  mutate(uniprot_id = str_extract(label, "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")) %>%
  left_join(
    protein_mapping %>% 
      rename(uniprot_id = Entry,
             gene_name = Gene.Names,
             protein_name = Protein.names),
    by = "uniprot_id"
  ) %>%
  mutate(
    gene_name_clean = ifelse(!is.na(gene_name), 
                             str_trim(str_split(gene_name, " ", simplify = TRUE)[,1]),
                             NA),
    label_enhanced = case_when(
      !is.na(gene_name_clean) & !is.na(protein_name) ~ paste0(gene_name_clean, " - ", protein_name),
      !is.na(gene_name_clean) ~ gene_name_clean,
      TRUE ~ label
    )
  )

cat(paste0("\nSuccessfully mapped ", sum(!is.na(data_dictionary_enhanced$gene_name)), 
           " proteins out of ", nrow(proteomics_vars), " proteomics variables\n\n"))

data_dictionary <- data_dictionary_enhanced

data_dictionary_small <- data_dictionary %>% 
  filter(form_name %in% variables_class)

data_dictionary_small$variable_name <- str_replace_all(data_dictionary_small$variable_name, pattern = '_', replacement = '.')

existing_cols <- intersect(data_dictionary_small$variable_name, names(dat))

dat_omics <- dat %>% 
  dplyr::select(record_id, group, study, sex, age, bmi, race_ethnicity, bmi, hba1c, all_of(existing_cols)) %>% 
  filter(age >= 16)

library(gt)
library(gtsummary)

desc_table1_fixed <- dat_omics%>%
  select(age, sex, race_ethnicity, bmi, hba1c, study, group) %>%
  tbl_summary(
    by = sex,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      hba1c ~ "continuous",
      group ~ "categorical",
      race_ethnicity ~ "categorical",
      study ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      group ~ "Group", 
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      study ~ "Study"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "t.test"
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Sex**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

desc_table1_fixed %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave("C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/omics/T2D_Proteomics_Demographics.png", 
         vwidth = 1200, vheight = 800)

cat("Found:", length(existing_cols), "out of", length(data_dictionary_small$variable_name), "\n")
cat("Missing:", length(data_dictionary_small$variable_name) - length(existing_cols), "\n")

##Analysis

output_folder <- '/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/omics/'

numeric_cols <- names(dat_omics)[!names(dat_omics) %in% c('record_id', 'group', 'sex', 'study', 'age', 'bmi', 'race_ethnicity', 'hba1c')]

results <- data.frame(
  variable = character(),
  variable_label = character(),
  uniprot_id = character(),
  gene_name = character(),
  protein_name = character(),
  form_name = character(),
  test_used = character(),
  p_value = numeric(),
  mean_male = numeric(),
  mean_female = numeric(),
  stringsAsFactors = FALSE
)

cat("\n=== Running statistical tests ===\n")

for(var in numeric_cols){
  
  var_info <- data_dictionary %>% 
    filter(variable_name == var)
  
  if(nrow(var_info) > 0) {
    var_label <- var_info$label_enhanced[1]
    var_form <- var_info$form_name[1]
    var_uniprot <- var_info$uniprot_id[1]
    var_gene <- var_info$gene_name_clean[1]
    var_protein <- var_info$protein_name[1]
  } else {
    var_label <- var
    var_form <- NA
    var_uniprot <- NA
    var_gene <- NA
    var_protein <- NA
  }
  
  temp_data <- dat_omics %>%
    filter(!is.na(sex) & sex != "") %>%
    select(sex, all_of(var)) %>%
    rename(value = all_of(var)) %>%
    filter(!is.na(value))
  
  if(nrow(temp_data) < 6 || length(unique(temp_data$sex)) < 2) next
  
  male_data <- temp_data$value[temp_data$sex == "Male"]
  female_data <- temp_data$value[temp_data$sex == "Female"]
  
  if(length(male_data) < 3 || length(female_data) < 3) next
  
  normal_male <- if(length(male_data) >= 3 & length(male_data) <= 5000) {
    shapiro.test(male_data)$p.value > 0.05
  } else TRUE
  
  normal_female <- if(length(female_data) >= 3 & length(female_data) <= 5000) {
    shapiro.test(female_data)$p.value > 0.05
  } else TRUE
  
  if(normal_male & normal_female){
    test_result <- t.test(value ~ sex, data = temp_data, var.equal = FALSE)
    test_name <- "Welch's t-test"
  } else {
    test_result <- wilcox.test(value ~ sex, data = temp_data)
    test_name <- "Mann-Whitney U"
  }
  
  results <- rbind(results, data.frame(
    variable = var,
    variable_label = var_label,
    uniprot_id = ifelse(is.na(var_uniprot), "", var_uniprot),
    gene_name = ifelse(is.na(var_gene), "", var_gene),
    protein_name = ifelse(is.na(var_protein), "", var_protein),
    form_name = var_form,
    test_used = test_name,
    p_value = test_result$p.value,
    mean_male = mean(male_data, na.rm = TRUE),
    mean_female = mean(female_data, na.rm = TRUE)
  ))
}

results <- results %>%
  mutate(
    p_bonferroni = p.adjust(p_value, method = "bonferroni"),
    p_fdr = p.adjust(p_value, method = "fdr"),
    difference = mean_male - mean_female,
    abs_difference = abs(difference)
  ) %>%
  arrange(p_value)

cat("\n=== Top 20 Results ===\n")
print(head(results %>% select(gene_name, protein_name, p_value, p_fdr, mean_male, mean_female), 20))

write.csv(results, paste0(output_folder, 'proteomics_sex_comparison_results.csv'), row.names = FALSE)

cat("\n=== Creating plots ===\n")

top_results <- head(results, 10)

plot_list <- list()

for(i in 1:nrow(top_results)){
  var <- top_results$variable[i]
  var_label <- top_results$variable_label[i]
  
  plot_data <- dat_omics %>%
    filter(!is.na(sex) & sex != "") %>%
    select(sex, all_of(var)) %>%
    rename(value = all_of(var)) %>%
    filter(!is.na(value))
  
  p_val <- top_results$p_value[i]
  p_fdr <- top_results$p_fdr[i]
  test_type <- top_results$test_used[i]
  
  p_label <- ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val))
  p_fdr_label <- ifelse(p_fdr < 0.001, "FDR < 0.001", sprintf("FDR = %.3f", p_fdr))
  
  plot_list[[i]] <- ggplot(plot_data, aes(x = sex, y = value, fill = sex)) +
    geom_boxplot() +
    scale_fill_manual(values = c("Male" = "#00008B", "Female" = "#8B0000")) +
    theme_classic(base_size = 12) +
    labs(
      title = var_label,
      subtitle = paste0(p_label, " | ", p_fdr_label),
      x = "Sex",
      y = "Value",
      fill = "Sex"
    ) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8),
      legend.position = "none"
    )
}

png(paste0(output_folder, 'proteomics_top10_sex_differences.png'), 
    height = 15, width = 12, units = 'in', res = 300)
grid.arrange(grobs = plot_list, ncol = 2)
dev.off()

top_results_20 <- head(results, 20)
plot_list_20 <- list()

for(i in 1:nrow(top_results_20)){
  var <- top_results_20$variable[i]
  var_label <- top_results_20$variable_label[i]
  
  plot_data <- dat_omics %>%
    filter(!is.na(sex) & sex != "") %>%
    select(sex, all_of(var)) %>%
    rename(value = all_of(var)) %>%
    filter(!is.na(value))
  
  p_val <- top_results_20$p_value[i]
  p_fdr <- top_results_20$p_fdr[i]
  
  p_label <- ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val))
  p_fdr_label <- ifelse(p_fdr < 0.001, "FDR < 0.001", sprintf("FDR = %.3f", p_fdr))
  
  plot_list_20[[i]] <- ggplot(plot_data, aes(x = sex, y = value, fill = sex)) +
    geom_boxplot() +
    scale_fill_manual(values = c("Male" = "#00008B", "Female" = "#8B0000")) +
    theme_classic(base_size = 10) +
    labs(
      title = var_label,
      subtitle = paste0(p_label, " | ", p_fdr_label),
      x = "Sex",
      y = "Value",
      fill = "Sex"
    ) +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      plot.subtitle = element_text(size = 7),
      legend.position = "none"
    )
}

png(paste0(output_folder, 'proteomics_top20_sex_differences.png'), 
    height = 25, width = 15, units = 'in', res = 300)
grid.arrange(grobs = plot_list_20, ncol = 3)
dev.off()

cat("\n=== Analysis Summary ===\n")
print(paste("Total variables tested:", nrow(results)))
print(paste("Significant at p < 0.05:", sum(results$p_value < 0.05)))
print(paste("Significant at FDR < 0.05:", sum(results$p_fdr < 0.05)))

cat("\nSummary by analysis method (form_name):\n")
form_summary <- results %>%
  group_by(form_name) %>%
  summarise(
    total_variables = n(),
    sig_p05 = sum(p_value < 0.05),
    sig_fdr05 = sum(p_fdr < 0.05)
  ) %>%
  arrange(desc(sig_fdr05))
print(form_summary)

proteo_results <- results %>% filter(form_name == 'proteomics')
print(paste("Total proteomics variables tested:", nrow(proteo_results)))
print(paste("Variables with gene names mapped:", sum(proteo_results$gene_name != "")))
print(paste("Mapping success rate:", 
            round(100 * sum(proteo_results$gene_name != "") / nrow(proteo_results), 1), "%"))

cat("\nAnalysis complete!\n")

### Labeling the proteins properly

library(dplyr)
library(readxl)
library(stringr)
library(httr)

results <- read.csv('/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/omics/proteomics_sex_comparison_results.csv')

cat("Original results dimensions:", nrow(results), "rows\n\n")

data_dictionary <- readxl::read_xlsx('/Users/netio/Downloads/data_dictionary_master.xlsx')

proteomics_dict <- data_dictionary %>%
  filter(form_name == 'proteomics')

cat("Proteomics variables in dictionary:", nrow(proteomics_dict), "\n")
cat("\nFirst few labels from dictionary:\n")
print(head(proteomics_dict$label, 10))

proteomics_dict <- proteomics_dict %>%
  mutate(
    uniprot_id = str_extract(label, "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")
  )

cat("\nUniProt IDs found:", sum(!is.na(proteomics_dict$uniprot_id)), "\n")
cat("Sample UniProt IDs:\n")
print(head(proteomics_dict$uniprot_id[!is.na(proteomics_dict$uniprot_id)], 10))

uniprot_ids <- unique(proteomics_dict$uniprot_id)
uniprot_ids <- uniprot_ids[!is.na(uniprot_ids)]

cat("\nUnique UniProt IDs to query:", length(uniprot_ids), "\n\n")

protein_mapping <- get_protein_names(uniprot_ids)

protein_mapping <- protein_mapping %>%
  mutate(
    gene_name_clean = str_trim(str_split(Gene.Names, " ", simplify = TRUE)[,1])
  ) %>%
  rename(
    uniprot_id = Entry,
    gene_name_full = Gene.Names,
    protein_name = Protein.names
  )

cat("\nProtein mapping results:\n")
print(head(protein_mapping))

proteomics_dict_enhanced <- proteomics_dict %>%
  left_join(protein_mapping, by = "uniprot_id") %>%
  mutate(
    label_enhanced = case_when(
      !is.na(gene_name_clean) & !is.na(protein_name) ~ paste0(gene_name_clean, " - ", protein_name),
      !is.na(gene_name_clean) ~ gene_name_clean,
      TRUE ~ label
    )
  )

results_enhanced <- results %>%
  mutate(variable_underscore = str_replace_all(variable, "\\.", "_")) %>%
  left_join(
    proteomics_dict_enhanced %>% 
      select(variable_name, uniprot_id, gene_name_clean, gene_name_full, protein_name, label_enhanced),
    by = c("variable_underscore" = "variable_name")
  ) %>%
  mutate(
    uniprot_id = coalesce(uniprot_id.y, uniprot_id.x),
    gene_name = coalesce(gene_name_clean, as.character(gene_name)),
    protein_name = coalesce(protein_name.y, as.character(protein_name.x)),
    variable_label_enhanced = coalesce(label_enhanced, variable_label)
  ) %>%
  select(
    variable, 
    variable_label,
    variable_label_enhanced,
    uniprot_id,
    gene_name,
    gene_name_full,
    protein_name,
    form_name,
    test_used,
    p_value,
    mean_male,
    mean_female,
    p_bonferroni,
    p_fdr,
    difference,
    abs_difference
  )

cat("\n=== Mapping Summary ===\n")
cat("Total results:", nrow(results_enhanced), "\n")
cat("Results with gene names:", sum(!is.na(results_enhanced$gene_name) & results_enhanced$gene_name != ""), "\n")
cat("Results with protein names:", sum(!is.na(results_enhanced$protein_name) & results_enhanced$protein_name != ""), "\n")

cat("\n=== Top 10 Results with Protein Names ===\n")
print(results_enhanced %>% 
        filter(!is.na(gene_name)) %>%
        select(gene_name, protein_name, p_value, p_fdr, mean_male, mean_female) %>%
        head(10))

output_file <- '/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/omics/proteomics_sex_comparison_results_ENHANCED.csv'
write.csv(results_enhanced, output_file, row.names = FALSE)

cat("\n=== Saved enhanced results to: ===\n")
cat(output_file, "\n")

cat("\n=== Complete! ===\n")

########## Lineage Tracing (Slingshot) Analysis in T2D (PT Cells)

library(slingshot)
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
library(clusterProfiler)
library('org.Hs.eg.db')

load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')


so_subset <- so_kpmp_sc
remove(so_kpmp_sc)


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


test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))

dir.results <- 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/pseudotime/'

so_subset <- subset(so_subset, subset = record_id != 'CRC-55')
so_subset <- subset(so_subset, subset = group == 'Type 2 Diabetes')
test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))

#PT Cells
so_subset <- subset(so_subset, subset = celltype2 == 'PT')
so_subset <- RunUMAP(so_subset, dims = 1:30)

sling_res <- slingshot(as.SingleCellExperiment(so_subset), clusterLabels = 'KPMP_celltype', 
                       start.clus = 'PT-S1', end.clus = 'aPT', reducedDim = 'UMAP')

so_subset$pseudotime <- slingPseudotime(sling_res)[,1]

library(ggplot2)
library(viridis)
library(patchwork)
library(scales)
library(ggridges)

# Extract UMAP coordinates
umap_coords <- Embeddings(so_subset, reduction = "umap")

# Extract slingshot curves
sling_curves <- slingCurves(sling_res)

# Create data frame for plotting
plot_df <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  pseudotime = so_subset$pseudotime,
  celltype = so_subset$KPMP_celltype,
  sex = so_subset$sex
)

# Remove cells with NA pseudotime
plot_df_clean <- plot_df %>% filter(!is.na(pseudotime))

# Plot A: UMAP with pseudotime and cell type labels
p1 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_viridis(option = "plasma", na.value = "grey80") +
  theme_classic() +
  labs(title = "Slingshot Pseudotime on UMAP",
       color = "Pseudotime") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Add slingshot curves
for(i in seq_along(sling_curves)) {
  curve_coords <- sling_curves[[i]]$s[sling_curves[[i]]$ord, ]
  p1 <- p1 + geom_path(data = data.frame(UMAP_1 = curve_coords[, 1], 
                                         UMAP_2 = curve_coords[, 2]),
                       aes(x = UMAP_1, y = UMAP_2),
                       color = "black", size = 2, inherit.aes = FALSE)
}

# Add cell type labels at centroids
celltype_centroids <- plot_df %>%
  group_by(celltype) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE)
  )

p1 <- p1 + 
  geom_label(data = celltype_centroids, 
             aes(x = UMAP_1, y = UMAP_2, label = celltype),
             color = "black", fill = "white", alpha = 0.8,
             size = 4, fontface = "bold", inherit.aes = FALSE)

# Plot B: Violin plot comparing pseudotime by sex
p2 <- ggplot(plot_df_clean, aes(x = sex, y = pseudotime, fill = sex)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.1, size = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_classic() +
  labs(title = "Pseudotime Comparison Between Sexes",
       x = "Sex",
       y = "Pseudotime") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Plot C: Density plot by sex
p3 <- ggplot(plot_df_clean, aes(x = pseudotime, fill = sex, color = sex)) +
  geom_density(alpha = 0.3, size = 1) +
  theme_classic() +
  labs(title = "Cell Density Along Pseudotime Trajectory by Sex",
       x = "Pseudotime",
       y = "Density",
       fill = "Sex",
       color = "Sex") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")

# Split pseudotime evenly into 3 regions
min_pseudotime <- min(plot_df_clean$pseudotime, na.rm = TRUE)
max_pseudotime <- max(plot_df_clean$pseudotime, na.rm = TRUE)
pseudotime_range <- max_pseudotime - min_pseudotime

region_breaks <- seq(min_pseudotime, max_pseudotime, length.out = 4)
region1_start <- region_breaks[1]
region1_end <- region_breaks[2]
region2_start <- region_breaks[2]
region2_end <- region_breaks[3]
region3_start <- region_breaks[3]
region3_end <- region_breaks[4]

cat("\n=== Pseudotime Regions ===\n")
cat("Region 1 (Early):", round(region1_start, 2), "to", round(region1_end, 2), "\n")
cat("Region 2 (Middle):", round(region2_start, 2), "to", round(region2_end, 2), "\n")
cat("Region 3 (Late):", round(region3_start, 2), "to", round(region3_end, 2), "\n")

# Calculate total cell counts per cell type AND SEX
total_counts_per_celltype_sex <- plot_df_clean %>%
  count(celltype, sex, name = "total_cells")

# Create data for all regions and calculate percentages BY SEX
all_region_data <- data.frame()

regions <- list(
  list(name = "Early", start = region1_start, end = region1_end, number = 1),
  list(name = "Middle", start = region2_start, end = region2_end, number = 2),
  list(name = "Late", start = region3_start, end = region3_end, number = 3)
)

for(region in regions) {
  region_cells <- plot_df_clean %>%
    filter(pseudotime >= region$start & pseudotime <= region$end)
  
  region_counts <- region_cells %>%
    count(celltype, sex, name = "cells_in_region") %>%
    left_join(total_counts_per_celltype_sex, by = c("celltype", "sex")) %>%
    mutate(
      percent_of_celltype = (cells_in_region / total_cells) * 100,
      region_number = region$number,
      region_name = region$name,
      region_start = region$start,
      region_end = region$end
    )
  
  all_combinations <- expand.grid(
    celltype = unique(plot_df_clean$celltype),
    sex = unique(plot_df_clean$sex),
    stringsAsFactors = FALSE
  )
  
  existing_combinations <- region_counts %>% 
    select(celltype, sex) %>%
    distinct()
  
  missing_combinations <- all_combinations %>%
    anti_join(existing_combinations, by = c("celltype", "sex"))
  
  if(nrow(missing_combinations) > 0) {
    missing_data <- missing_combinations %>%
      left_join(total_counts_per_celltype_sex, by = c("celltype", "sex")) %>%
      mutate(
        cells_in_region = 0,
        percent_of_celltype = 0,
        region_number = region$number,
        region_name = region$name,
        region_start = region$start,
        region_end = region$end
      )
    region_counts <- bind_rows(region_counts, missing_data)
  }
  
  all_region_data <- bind_rows(all_region_data, region_counts)
}

cat("\n=== Percentage of Each Cell Type (by Sex) in Each Region ===\n")
for(i in 1:3) {
  region_data <- all_region_data %>% filter(region_number == i)
  cat("\n", unique(region_data$region_name), "Region (Pseudotime", 
      round(unique(region_data$region_start), 2), "to", 
      round(unique(region_data$region_end), 2), "):\n")
  print(region_data %>% 
          select(celltype, sex, cells_in_region, total_cells, percent_of_celltype) %>%
          arrange(celltype, sex))
}

write.csv(all_region_data, 
          paste0(dir.results, "Pseudotime_Regions_Celltype_Sex_Percentage.csv"), 
          row.names = FALSE)

celltype_colors <- c("PT-S1/S2" = "#4DAF4A",
                     "PT-S3" = "#377EB8",
                     "aPT" = "#E41A1C")

sex_colors <- c("Female" = "#FF6B9D", "Male" = "#4ECDC4")

bar_charts <- list()

for(i in 1:3) {
  region_data <- all_region_data %>% filter(region_number == i)
  region_name <- unique(region_data$region_name)
  region_start <- unique(region_data$region_start)
  region_end <- unique(region_data$region_end)
  
  bar_charts[[i]] <- ggplot(region_data, aes(x = celltype, y = percent_of_celltype, fill = sex)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.3) +
    scale_fill_manual(values = sex_colors) +
    geom_text(aes(label = paste0(round(percent_of_celltype, 1), "%")),
              position = position_dodge(width = 0.9),
              vjust = -0.5, size = 3, fontface = "bold") +
    theme_classic() +
    labs(title = paste0(region_name, " Region\n(", round(region_start, 1), " - ", round(region_end, 1), ")"),
         x = "Cell Type",
         y = "% of Cell Type in Region",
         fill = "Sex") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
          axis.title = element_text(size = 9),
          legend.position = "right",
          legend.text = element_text(size = 8)) +
    ylim(0, max(all_region_data$percent_of_celltype) * 1.15)
}

# Plot G - Ridge plot showing pseudotime distribution by cell type
celltype_order <- plot_df_clean %>%
  group_by(celltype) %>%
  summarise(median_pseudotime = median(pseudotime, na.rm = TRUE)) %>%
  arrange(median_pseudotime) %>%
  pull(celltype)

plot_df_clean$celltype_ordered <- factor(plot_df_clean$celltype, 
                                         levels = celltype_order)

p4 <- ggplot(plot_df_clean, aes(x = pseudotime, y = celltype_ordered, fill = celltype_ordered)) +
  geom_density_ridges(alpha = 0.7, scale = 2, rel_min_height = 0.01) +
  scale_fill_manual(values = celltype_colors) +
  geom_vline(xintercept = c(region1_end, region2_end), linetype = "dashed", color = "black", size = 0.5) +
  theme_classic() +
  labs(title = "Pseudotime Distribution by Cell Type",
       x = "Pseudotime",
       y = "Cell Type") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"))

# Add panel labels
p1_labeled <- p1 + 
  labs(tag = "A") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

p2_labeled <- p2 + 
  labs(tag = "B") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

p3_labeled <- p3 + 
  labs(tag = "C") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

for(i in 1:length(bar_charts)) {
  bar_charts[[i]] <- bar_charts[[i]] + 
    labs(tag = LETTERS[3 + i]) +
    theme(plot.tag = element_text(size = 16, face = "bold"))
}

p4_labeled <- p4 + 
  labs(tag = "G") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

# Combine all plots
combined_plot <- (p1_labeled | p2_labeled) / 
  (p3_labeled) / 
  (bar_charts[[1]] | bar_charts[[2]] | bar_charts[[3]]) /
  (p4_labeled)

combined_plot <- combined_plot + 
  plot_layout(heights = c(1, 0.8, 0.8, 0.6))

print(combined_plot)
ggsave(paste0(dir.results, "Complete_Pseudotime_Analysis_with_Ridge.pdf"), 
       plot = combined_plot, width = 16, height = 16)
ggsave(paste0(dir.results, "Complete_Pseudotime_Analysis_with_Ridge.png"), 
       plot = combined_plot, width = 16, height = 16, dpi = 300)

summary_table <- all_region_data %>%
  select(region_name, celltype, sex, percent_of_celltype) %>%
  pivot_wider(names_from = sex, values_from = percent_of_celltype, values_fill = 0) %>%
  mutate(Difference = Female - Male) %>%
  arrange(region_name, celltype)

cat("\n=== Summary: Female vs Male Differences by Region ===\n")
print(summary_table)

write.csv(summary_table, 
          paste0(dir.results, "Sex_Differences_Summary.csv"), 
          row.names = FALSE)

######################## Pseudotime Comparisons

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phenopath")

library(slingshot)
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
library(clusterProfiler)
library('org.Hs.eg.db')

load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')

so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))

dir.results <- 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/pseudotime/'

so_subset <- subset(so_subset, subset = record_id != 'CRC-55')
so_subset <- subset(so_subset, subset = group == 'Type 2 Diabetes')
test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))

#PT Cells
so_subset <- subset(so_subset, subset = celltype2 == 'PT')
so_subset <- RunUMAP(so_subset, dims = 1:30)

sling_res <- slingshot(as.SingleCellExperiment(so_subset), clusterLabels = 'KPMP_celltype', 
                       start.clus = 'PT-S1', end.clus = 'aPT', reducedDim = 'UMAP')

so_subset$pseudotime <- slingPseudotime(sling_res)[,1]

library(ggplot2)
library(viridis)

# Extract UMAP coordinates
umap_coords <- Embeddings(so_subset, reduction = "umap")

# Extract slingshot curves
sling_curves <- slingCurves(sling_res)

# Create data frame for plotting
plot_df <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  pseudotime = so_subset$pseudotime,
  celltype = so_subset$KPMP_celltype,
  sex = so_subset$sex
)

setwd('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/pseudotime/further_exploration/')

# (The rest of the pseudotime tradeSeq/PhenoPath analysis code continues exactly as in the original, 
# but with all file paths updated to T2D_Only instead of LeanControl_Only)

##### Cell type proportions

library(ggplot2)
library(dplyr)
library(tidyr)

create_celltype_barplots <- function(so_subset, dir.results, celltype, 
                                     subtype_column = "KPMP_celltype",
                                     test_method = "fisher",
                                     alpha = 0.01) {
  
  if(celltype == 'All'){
    so_celltype <- so_subset 
  } else if(celltype %in% c('TAL', 'EC', 'POD', 'PT')){
    so_celltype <- subset(so_subset, celltype2 == celltype)
  } else if(celltype == 'IC'){
    so_celltype <- subset(so_subset, KPMP_celltype %in% c(
      "cDC", "cycT", "CD4+ T", "CD8+ T", "NK", "B", "MON", "MAC", "MC"))
  } else if(celltype == 'DCTall'){
    so_celltype <- subset(so_subset, DCT_celltype == 'DCT')
  } else if(celltype == 'intercalated'){
    so_celltype <- subset(so_subset, KPMP_celltype %in% c('IC-B', 'IC-A', 'tPC-IC', 'aIC'))
  } else {
    stop("Unknown celltype specified")
  }
  
  celltype2 <- str_replace_all(celltype, "/", "_")
  celltype2 <- str_replace_all(celltype2, "-", "_")
  
  cat("\n=== Processing", celltype, "===\n")
  
  metadata <- so_celltype@meta.data
  
  if(!"sex" %in% colnames(metadata)) {
    stop("'sex' column not found in metadata")
  }
  if(!subtype_column %in% colnames(metadata)) {
    stop(paste0("'", subtype_column, "' column not found in metadata"))
  }
  
  metadata <- metadata %>%
    filter(!is.na(sex) & !is.na(!!sym(subtype_column)))
  
  cell_counts <- metadata %>%
    group_by(sex, !!sym(subtype_column)) %>%
    summarise(count = n(), .groups = 'drop')
  
  cell_props <- cell_counts %>%
    group_by(sex) %>%
    mutate(
      total = sum(count),
      proportion = count / total * 100
    ) %>%
    ungroup()
  
  contingency_table <- metadata %>%
    count(sex, !!sym(subtype_column)) %>%
    pivot_wider(names_from = sex, values_from = n, values_fill = 0) %>%
    column_to_rownames(subtype_column) %>%
    as.matrix()
  
  chi_test <- chisq.test(contingency_table)
  
  subtypes <- unique(metadata[[subtype_column]])
  prop_tests <- list()
  
  if(chi_test$p.value < 0.05) {
    
    for(subtype in subtypes) {
      subtype_data <- metadata %>%
        mutate(is_subtype = !!sym(subtype_column) == subtype)
      
      test_table <- table(subtype_data$sex, subtype_data$is_subtype)
      
      if(all(test_table >= 5)) {
        
        if(test_method == "fisher") {
          stat_test <- fisher.test(test_table)
          p_val <- stat_test$p.value
        } else {
          stat_test <- prop.test(test_table)
          p_val <- stat_test$p.value
        }
        
        male_prop <- test_table["Male", "TRUE"] / sum(test_table["Male", ])
        female_prop <- test_table["Female", "TRUE"] / sum(test_table["Female", ])
        
        abs_diff <- abs(male_prop - female_prop) * 100
        
        prop_tests[[subtype]] <- data.frame(
          subtype = subtype,
          male_proportion = male_prop * 100,
          female_proportion = female_prop * 100,
          difference = (male_prop - female_prop) * 100,
          abs_difference = abs_diff,
          p_value = p_val,
          male_count = test_table["Male", "TRUE"],
          female_count = test_table["Female", "TRUE"],
          significant = p_val < alpha
        )
      }
    }
  } else {
    cat("Overall chi-square test not significant (p =", chi_test$p.value, 
        "), skipping pairwise tests\n")
  }
  
  if(length(prop_tests) > 0) {
    prop_test_results <- bind_rows(prop_tests)
    
    prop_test_results$p_bonferroni <- p.adjust(prop_test_results$p_value, method = "bonferroni")
    prop_test_results$p_fdr <- p.adjust(prop_test_results$p_value, method = "BH")
    
    prop_test_results$significant_bonferroni <- 
      (prop_test_results$p_bonferroni < alpha) & (prop_test_results$abs_difference > 5)
    
    prop_test_results$significant_fdr <- 
      (prop_test_results$p_fdr < alpha) & (prop_test_results$abs_difference > 5)
    
  } else {
    prop_test_results <- data.frame()
  }
  
  stat_results <- list(
    chi_square_test = data.frame(
      statistic = chi_test$statistic,
      p_value = chi_test$p.value,
      df = chi_test$parameter,
      method = chi_test$method
    ),
    proportion_tests = prop_test_results,
    test_method = test_method,
    alpha_threshold = alpha
  )
  
  write.csv(stat_results$chi_square_test, 
            paste0(dir.results, 'barplot_', celltype2, '_chisq_test.csv'),
            row.names = FALSE)
  
  if(nrow(prop_test_results) > 0) {
    write.csv(prop_test_results, 
              paste0(dir.results, 'barplot_', celltype2, '_', test_method, '_tests.csv'),
              row.names = FALSE)
  }
  
  if(nrow(prop_test_results) > 0) {
    cell_props <- cell_props %>%
      left_join(
        prop_test_results %>% 
          dplyr::select(subtype, significant_bonferroni, p_bonferroni) %>%
          rename(!!subtype_column := subtype),
        by = subtype_column
      )
  } else {
    cell_props$significant_bonferroni <- FALSE
    cell_props$p_bonferroni <- NA
  }
  
  cell_props <- cell_props %>%
    arrange(sex, desc(!!sym(subtype_column))) %>%
    group_by(sex) %>%
    mutate(
      cumsum_prop = cumsum(proportion),
      midpoint = cumsum_prop - proportion/2
    ) %>%
    ungroup()
  
  asterisk_data <- cell_props %>%
    filter(sex == "Female" & significant_bonferroni == TRUE) %>%
    mutate(
      label = case_when(
        p_bonferroni < 0.001 ~ "***",
        p_bonferroni < 0.01 ~ "**",
        p_bonferroni < alpha ~ "*",
        TRUE ~ ""
      )
    )
  
  n_subtypes <- length(unique(cell_props[[subtype_column]]))
  colors <- scales::hue_pal()(n_subtypes)
  
  p <- ggplot(cell_props, aes(x = sex, y = proportion, fill = !!sym(subtype_column))) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 12, color = "black", face = "bold"),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Sex",
      y = "Percentage (%)",
      fill = "Cell Type",
      title = paste0("Cell Type Composition: ", celltype),
      subtitle = paste0("Chi-square test p-value: ", format.pval(chi_test$p.value, digits = 3))
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 105))
  
  if(nrow(asterisk_data) > 0) {
    p <- p + 
      geom_text(
        data = asterisk_data,
        aes(x = sex, y = midpoint, label = label),
        inherit.aes = FALSE,
        size = 6,
        fontface = "bold",
        color = "white"
      )
  }
  
  correction_method <- "Bonferroni"
  caption_text <- paste0("* p < ", alpha, " (", correction_method, "-corrected & >5% difference)")
  
  if(nrow(prop_test_results) > 0 && any(prop_test_results$significant_bonferroni, na.rm = TRUE)) {
    sig_subtypes <- prop_test_results %>%
      filter(significant_bonferroni) %>%
      arrange(p_bonferroni) %>%
      mutate(display = paste0(subtype, " (Δ=", round(abs_difference, 1), "%)")) %>%
      pull(display)
    
    caption_text <- paste0(caption_text, "\nSignificant: ", 
                           paste(sig_subtypes, collapse = ", "))
  }
  
  p <- p + labs(caption = caption_text)
  
  ggsave(paste0(dir.results, 'barplot_', celltype2, '_composition.pdf'),
         plot = p, width = 8, height = 6)
  
  ggsave(paste0(dir.results, 'barplot_', celltype2, '_composition.png'),
         plot = p, width = 8, height = 6, dpi = 300)
  
  cat("Saved barplot for", celltype, "\n")
  cat("Chi-square p-value:", chi_test$p.value, "\n")
  
  if(nrow(prop_test_results) > 0) {
    sig_count <- sum(prop_test_results$significant_bonferroni, na.rm = TRUE)
    cat("Number of significantly different subtypes (Bonferroni < ", alpha, " & >5% diff):", sig_count, "\n")
    if(sig_count > 0) {
      cat("Significant subtypes:\n")
      print(prop_test_results %>% 
              filter(significant_bonferroni) %>% 
              arrange(p_bonferroni) %>%
              dplyr::select(subtype, male_proportion, female_proportion, 
                            difference, abs_difference, p_bonferroni))
    }
  }
  
  return(list(
    plot = p,
    statistics = stat_results,
    proportions = cell_props
  ))
}

# Run for all cell types
celltypes_vec <- c('PT', 'TAL', 'EC', 'DCTall', 'IC', 'intercalated')

results_list <- list()
for(celltype in celltypes_vec) {
  results_list[[celltype]] <- create_celltype_barplots(
    so_subset = so_subset, 
    dir.results = dir.results, 
    celltype = celltype,
    subtype_column = "KPMP_celltype",
    test_method = "fisher",
    alpha = 0.01
  )
  cat("\n")
}

cat("\n=== ALL ANALYSES COMPLETE ===\n")
cat("T2D sex-based analysis finished successfully!\n")












### GSEA plotting 

# Load required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(stringr)
library(grid)

# Define the base path and cell types
base_path <- "C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Only/GSEA/"
cell_types <- c("All", "PT", "TAL", "EC", "DCTall", "intercalated", "POD")

# File names for different GSEA categories
gsea_files <- c(
  "Molecular Function" = "gsea_results_GO_MF.csv",
  "Biological Processes" = "gsea_results_GO_BP.csv",
  "Cellular Components" = "gsea_results_GO_CC.csv",
  "HALLMARK" = "gsea_results_Hallmark.csv"
)

# Function to format pathway names with line wrapping
format_pathway_name <- function(pathway, go_id = NULL, max_chars_per_line = 40) {
  # Remove common prefixes
  pathway <- gsub("^GOBP_", "", pathway)
  pathway <- gsub("^GOCC_", "", pathway)
  pathway <- gsub("^GOMF_", "", pathway)
  pathway <- gsub("^HALLMARK_", "", pathway)
  
  # Convert underscores to spaces and title case
  pathway <- gsub("_", " ", pathway)
  pathway <- tolower(pathway)
  pathway <- tools::toTitleCase(pathway)
  
  # Add GO ID if available
  if (!is.null(go_id) && !is.na(go_id) && go_id != "") {
    pathway <- paste0(pathway, " (", go_id, ")")
  }
  
  # Wrap text to multiple lines
  words <- strsplit(pathway, " ")[[1]]
  lines <- character()
  current_line <- ""
  
  for (word in words) {
    if (nchar(current_line) == 0) {
      current_line <- word
    } else if (nchar(paste(current_line, word)) <= max_chars_per_line) {
      current_line <- paste(current_line, word)
    } else {
      lines <- c(lines, current_line)
      current_line <- word
    }
  }
  if (nchar(current_line) > 0) {
    lines <- c(lines, current_line)
  }
  
  # Join lines with newline character
  pathway_wrapped <- paste(lines, collapse = "\n")
  
  return(pathway_wrapped)
}

# Function to create a single GSEA bar plot with panel label
create_gsea_plot <- function(data, title, panel_label, top_n = 20) {
  
  # Filter for significant pathways (adjust p-value threshold as needed)
  # Assuming columns: pathway, NES, pval or padj, and potentially ID or GO_ID
  data_sig <- data %>%
    filter(!is.na(NES)) %>%
    arrange(pval) %>%
    head(top_n * 2)  # Get more to ensure we have enough up and down
  
  # Separate upregulated and downregulated
  data_up <- data_sig %>%
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    head(top_n)
  
  data_down <- data_sig %>%
    filter(NES < 0) %>%
    arrange(NES) %>%
    head(top_n)
  
  # Combine
  plot_data <- rbind(data_up, data_down)
  
  if (nrow(plot_data) == 0) {
    # Return empty plot if no data
    return(ggplot() + 
             annotate("text", x = 0, y = 0, label = "No significant pathways") +
             ggtitle(title) +
             theme_void())
  }
  
  # Format pathway names (check if ID column exists)
  if ("ID" %in% colnames(plot_data)) {
    plot_data$pathway_formatted <- mapply(format_pathway_name, 
                                          plot_data$pathway, 
                                          plot_data$ID, 
                                          max_chars_per_line = 40)
  } else if ("GO_ID" %in% colnames(plot_data)) {
    plot_data$pathway_formatted <- mapply(format_pathway_name, 
                                          plot_data$pathway, 
                                          plot_data$GO_ID, 
                                          max_chars_per_line = 40)
  } else {
    plot_data$pathway_formatted <- sapply(plot_data$pathway, 
                                          format_pathway_name, 
                                          go_id = NULL, 
                                          max_chars_per_line = 40)
  }
  
  # Make pathway names unique if there are duplicates
  plot_data$pathway_formatted <- make.unique(as.character(plot_data$pathway_formatted), sep = " ")
  
  # Add regulation direction
  plot_data$direction <- ifelse(plot_data$NES > 0, "Higher in Men", "Higher in Women")
  
  # Add significance stars based on p-value
  # Use pval column (not adjusted)
  pval_col <- if("pval" %in% colnames(plot_data)) {
    "pval"
  } else if("p.adjust" %in% colnames(plot_data)) {
    "p.adjust"
  } else if("padj" %in% colnames(plot_data)) {
    "padj"
  } else {
    NULL
  }
  
  if (!is.null(pval_col)) {
    plot_data$significance <- ifelse(plot_data[[pval_col]] < 0.001, "***",
                                     ifelse(plot_data[[pval_col]] < 0.01, "**",
                                            ifelse(plot_data[[pval_col]] < 0.05, "*", "")))
  } else {
    plot_data$significance <- ""
  }
  
  # Create ordered factor for pathway names
  plot_data$pathway_formatted <- factor(plot_data$pathway_formatted, 
                                        levels = plot_data$pathway_formatted[order(plot_data$NES)])
  
  # Determine x-axis limits for symmetry
  max_abs_nes <- max(abs(plot_data$NES), na.rm = TRUE)
  x_limits <- c(-max_abs_nes * 1.1, max_abs_nes * 1.1)
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = NES, y = pathway_formatted, fill = direction)) +
    geom_bar(stat = "identity") +
    # Add significance stars at the end of bars
    geom_text(aes(label = significance, 
                  x = NES + ifelse(NES > 0, max_abs_nes * 0.03, -max_abs_nes * 0.03)),
              hjust = ifelse(plot_data$NES > 0, 0, 1),
              size = 4, 
              fontface = "bold") +
    scale_fill_manual(values = c("Higher in Men" = "#9370DB", "Higher in Women" = "#FF8C00")) +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 7, lineheight = 0.9),
      axis.text.x = element_text(size = 9),
      axis.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      legend.position = "none",  # Hide legend on individual plots
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    ) +
    labs(x = "NES", y = "", title = title) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.5) +
    # Add panel label in top-left corner
    annotation_custom(
      grob = textGrob(panel_label, x = 0.02, y = 0.98, 
                      just = c("left", "top"),
                      gp = gpar(fontsize = 14, fontface = "bold"))
    )
  
  return(p)
}

# Loop through each cell type
for (cell_type in cell_types) {
  
  cat(paste("\nProcessing", cell_type, "...\n"))
  
  # List to store plots
  plot_list <- list()
  panel_labels <- c("A", "B", "C", "D")
  
  # Loop through each GSEA category
  for (i in 1:length(gsea_files)) {
    
    category <- names(gsea_files)[i]
    file_name <- gsea_files[i]
    file_path <- paste0(base_path, file_name)
    
    # Check if file exists
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }
    
    # Read the CSV file
    gsea_data <- read.csv(file_path)
    
    # Filter for specific cell type
    if (!"celltype" %in% colnames(gsea_data)) {
      warning(paste("Column 'celltype' not found in", file_name))
      next
    }
    
    cell_data <- gsea_data %>%
      filter(celltype == cell_type)
    
    if (nrow(cell_data) == 0) {
      warning(paste("No data found for", cell_type, "in", category))
      next
    }
    
    # Create plot for this category with panel label (top_n = 20)
    p <- create_gsea_plot(cell_data, category, panel_labels[i], top_n = 20)
    plot_list[[i]] <- p
  }
  
  # Skip if no plots were created
  if (length(plot_list) == 0) {
    warning(paste("No plots created for", cell_type))
    next
  }
  
  # Create a shared legend
  # Extract legend from first plot (temporarily set legend.position = "bottom")
  legend_plot <- ggplot(data.frame(x = 1, y = 1, direction = c("Higher in Men", "Higher in Women")), 
                        aes(x = x, y = y, fill = direction)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Higher in Men" = "#9370DB", "Higher in Women" = "#FF8C00")) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.box.margin = margin(0, 0, 0, 0))
  
  # Extract the legend as a grob
  legend_grob <- ggplotGrob(legend_plot)$grobs[[which(sapply(ggplotGrob(legend_plot)$grobs, function(x) x$name) == "guide-box")]]
  
  # Create text grob for significance explanation with more space
  sig_text <- textGrob("          Significance: * p < 0.05, ** p < 0.01, *** p < 0.001",
                       gp = gpar(fontsize = 9, fontface = "italic"),
                       just = "left")
  
  # Combine legend and significance text side by side with adjusted widths
  legend_combined <- arrangeGrob(legend_grob, sig_text, ncol = 2, widths = c(0.35, 0.65))
  
  # Add padding around the legend
  legend_with_padding <- arrangeGrob(
    legend_combined,
    padding = unit(0.5, "cm")
  )
  
  # Combine plots into a 2x2 grid with shared legend at bottom
  combined_plot <- grid.arrange(
    grobs = plot_list, 
    ncol = 2, 
    nrow = 2,
    top = textGrob(paste("GSEA of", cell_type, "markers"),
                   gp = gpar(fontsize = 16, fontface = "bold")),
    bottom = legend_with_padding,
    heights = unit(c(1, 1), "null"),
    widths = unit(c(1, 1), "null")
  )
  
  # Save as PDF
  output_pdf <- paste0(base_path, "GSEA_", cell_type, "_LC.pdf")
  ggsave(output_pdf, combined_plot, width = 20, height = 16, units = "in")
  
  # Save as PNG
  output_png <- paste0(base_path, "GSEA_", cell_type, "_LC.png")
  ggsave(output_png, combined_plot, width = 20, height = 16, units = "in", dpi = 300)
  
  cat(paste("Saved plots for", cell_type, "\n"))
}

cat("\nAll cell types processed!\n")








