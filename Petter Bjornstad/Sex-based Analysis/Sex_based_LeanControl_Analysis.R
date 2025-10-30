#### Sex-Specific Analyses in Lean Controls




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







load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')

so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

#dat_groups <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_GroupAssignments.txt')
#dat_groups <- dat_groups %>% filter(group2 %in% c('Lean Control', 'T2D-No SGLTi2'))

#so_subset <- subset(so_subset, record_id == dat_groups$record_id)
test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))








dir.results <- 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/'


so_subset <- subset(so_subset, subset = record_id != 'CRC-55')
so_subset <- subset(so_subset, subset = group == 'Lean_Control')
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
    # Skip categorical p-values if they cause issues
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Group**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

# Save version with epic
desc_table1_fixed %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave(paste0(dir.results, "LC_demographics.png"), 
         vwidth = 1200, vheight = 800)












#### Celltype Analyses



## Barplots and UMAPs 

dir.results <- '/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/'



png(paste0(dir.results, 'LC_UMAP.png'), width = 1200, height = 1200)
DimPlot(so_subset) + 
  xlab("UMAP 1") + 
  ylab("UMAP 2") +
  theme(
    axis.title = element_text(size = 18),  # Increase axis label size
    axis.text = element_text(size = 16),   # Increase axis tick label size
    legend.text = element_text(size = 13), # Increase legend text size
    legend.title = element_text(size = 15) # Increase legend title size
  )
dev.off()






# Extract cell type and group information
cell_data <- data.frame(
  cell_type = so_subset$celltype2,
  group = so_subset$sex  # or whatever your grouping variable is
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
    x = "Sex in Lean Controls",
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
# Calculate position above each bar for each cell type
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
    x = "Sex in Lean Controls",
    y = "Proportion (%)",
    fill = "Cell Type"
  ) +
  ylim(0, 105) +  # Extend y-axis to fit asterisks
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

ggsave(paste0(dir.results, 'cell_type_proportions_asterisks_above.png'), 
       plot = p2, width = 8, height = 6, dpi = 300)

# Print the statistical results
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
      x = "Sex in Lean Controls",
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
  # Calculate position above each bar for each cell type
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
      x = "Sex in Lean Controls",
      y = "Proportion (%)",
      fill = "Cell Type"
    ) +
    ylim(0, 105) +  # Extend y-axis to fit asterisks
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12)
    )
  
  ggsave(paste0(dir.results, celltype, '_proportions_asterisks_above.png'), 
         plot = p2, width = 8, height = 6, dpi = 300)
  
  # Print the statistical results
  print(stat_results)
  
  
  print(celltype)
  
  
  
  
  
  
  
  
  
  
  
  
  
}










##################aPT analysis


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
# Adjust sex level name as needed

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

# Print model summary
print(summary(model))

# Save
ggsave(paste0(dir.results, "aPT_dTAL_interaction_plot.png"), plot = p, 
       width = 8, height = 6, dpi = 300)




pathway_data <- data.table::fread("/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/omics/metabolomics_pathways/pathway_results.csv")


# Filter for significant pathways
significant_pathways <- pathway_data[pathway_data$`Raw p` < 0.05, ]

# Create shortened labels
shorten_label <- function(label, max_length = 30) {
  ifelse(nchar(label) > max_length, 
         paste0(substr(label, 1, max_length), "..."), 
         label)
}

significant_pathways$label_short <- shorten_label(significant_pathways$V1)

png('/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/omics/metabolomics_pathways/metabolomics_pathways.png', 
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
LC_NEBULA_Analysis <- function(so_subset, dir.results, celltype, genes){
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
  #FOR LOOP OVER VARIABLES 
 
    counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round
    count_gene <- counts_path
    meta_gene <- subset(so_celltype)@meta.data
    
    
    
    complete_idx <- complete.cases(meta_gene$sex)
    cat("Cells with complete data:", sum(complete_idx), "\n")
    
    # Step 3: Filter all your data to only include cells with complete predictor data
    meta_gene <- meta_gene[complete_idx, ]
    count_gene <- count_gene[, complete_idx]  # Note: subsetting columns for cells
    
    
    num_cells <- nrow(meta_gene)
    num_part <- unique(meta_gene$record_id) %>% length()
    
    tmp_df <- meta_gene %>% dplyr::select(record_id, group, sex) %>% filter(!duplicated(record_id))
    
    num_male <- tmp_df %>% filter(sex == 'Male') %>% nrow()
    num_female <- tmp_df %>% filter(sex == 'Female') %>% nrow()
    
    # Step 4: Create prediction matrix from the complete data
    pred_gene <- model.matrix(~sex, data = meta_gene)
    
    
    # library <- meta_gene$library_size
    library <- meta_gene$pooled_offset
    data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)
    
    if (is.null(data_g_gene)) {
      data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
    }
    
    #With offset
    result <- nebula(count = data_g_gene$count, id = data_g_gene$id, 
                     pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)
    
    
    
    #Make dataframe of final results
    full_results <- as.data.frame(result)
    #Calculate number of genes filtered out for low expression 
    
    
    full_results$num_cells <- num_cells
    full_results$num_male <- num_male
    full_results$num_female <- num_female
    
    write.table(full_results,paste0(dir.results,"NEBULA_", 
                                    celltype2, "_cells__LC_pooledoffset.csv"),
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
  LC_NEBULA_Analysis(so_subset = so_subset, dir.results = dir.results, celltype = celltypes_vec[i], genes = gene_list)
  print(celltypes_vec[i])
}







### Full Analysis (all Genes)



#function
LC_NEBULA_Analysis_full <- function(so_subset, dir.results, celltype, genes){
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
  
  
  
  
 # full_analysis <- FindVariableFeatures(so_celltype, selection.method = "vst", nfeatures = 2000)
#  hvgs <- VariableFeatures(full_analysis)
#  so_celltype <- subset(so_celltype, features = hvgs)
  
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  
  
  print(paste0(celltype2, ' is running.'))
  #FOR LOOP OVER VARIABLES 
  
  counts_path <- round(GetAssayData(so_celltype, layer = "counts")) # load counts and round
  count_gene <- counts_path
  meta_gene <- subset(so_celltype)@meta.data
  
  
  
  complete_idx <- complete.cases(meta_gene$sex)
  cat("Cells with complete data:", sum(complete_idx), "\n")
  
  # Step 3: Filter all your data to only include cells with complete predictor data
  meta_gene <- meta_gene[complete_idx, ]
  count_gene <- count_gene[, complete_idx]  # Note: subsetting columns for cells
  
  
  num_cells <- nrow(meta_gene)
  num_part <- unique(meta_gene$record_id) %>% length()
  
  tmp_df <- meta_gene %>% dplyr::select(record_id, group, sex) %>% filter(!duplicated(record_id))
  
  num_male <- tmp_df %>% filter(sex == 'Male') %>% nrow()
  num_female <- tmp_df %>% filter(sex == 'Female') %>% nrow()
  
  # Step 4: Create prediction matrix from the complete data
  pred_gene <- model.matrix(~sex, data = meta_gene)
  
  
  # library <- meta_gene$library_size
  library <- meta_gene$pooled_offset
  data_g_gene <- group_cell(count = count_gene, id = meta_gene$kit_id, pred = pred_gene,offset=library)
  
  if (is.null(data_g_gene)) {
    data_g_gene <- list(count = count_gene, id = meta_gene$kit_id, pred = pred_gene, offset = library)
  }
  
  #With offset
  result <- nebula(count = data_g_gene$count, id = data_g_gene$id, 
                   pred = data_g_gene$pred, ncore = 1, reml=T,model="NBLMM",output_re = T,covariance=T,offset=data_g_gene$library)
  
  
  
  #Make dataframe of final results
  full_results <- as.data.frame(result)
  #Calculate number of genes filtered out for low expression 
  
  
  full_results$num_cells <- num_cells
  full_results$num_male <- num_male
  full_results$num_female <- num_female
  
  write.table(full_results,paste0(dir.results,"Full_NEBULA_", 
                                  celltype2, "_cells__LC_pooledoffset.csv"),
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
  LC_NEBULA_Analysis_full(so_subset = so_subset, dir.results = dir.results, celltype = celltypes_vec[i], genes = gene_list)
  print(celltypes_vec[i])
}









########### Volcano plots 

library(ggplot2)
library(ggbreak)
library(dplyr)

variable_names <- c('All', 'PT', 'TAL', 'EC', 'IC', 'POD', 'DCTall', 'intercalated')

for(i in c(1:length(variable_names))){
  
  sig_markers <- data.table::fread(paste0('/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/Full_NEBULA_', 
                                          variable_names[i], '_cells__LC_pooledoffset.csv'))
  
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
  
  # Automatically detect if axis break is needed
  logfc_range <- range(tmp_df$LogFC, na.rm = TRUE)
  logfc_q95 <- quantile(abs(tmp_df$LogFC), 0.95, na.rm = TRUE)
  logfc_max <- max(abs(tmp_df$LogFC), na.rm = TRUE)
  
  # Determine if we need a break (outliers beyond 95th percentile create large gaps)
  needs_break <- logfc_max > (logfc_q95 * 2)
  
  # Calculate break points if needed
  if(needs_break){
    # Find gaps in the data
    sorted_logfc <- sort(abs(tmp_df$LogFC))
    # Calculate gaps between consecutive values
    gaps <- diff(sorted_logfc)
    # Find the largest gap in the outer 10% of data
    gap_threshold <- quantile(abs(tmp_df$LogFC), 0.90, na.rm = TRUE)
    outer_data_idx <- which(sorted_logfc > gap_threshold)
    
    if(length(outer_data_idx) > 1){
      outer_gaps <- gaps[outer_data_idx[-length(outer_data_idx)]]
      max_gap_idx <- which.max(outer_gaps)
      
      # Set break points around the largest gap
      break_start <- sorted_logfc[outer_data_idx[max_gap_idx]] + 0.1
      break_end <- sorted_logfc[outer_data_idx[max_gap_idx + 1]] - 0.1
      
      # Ensure breaks are reasonable
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
  
  # Making graph
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
    
    # Add break if needed
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
    
    # Add break if needed
    if(use_break){
      tmp_graph <- tmp_graph + scale_x_break(breaks = c(break_start, break_end), scales = 0.3)
    }
  }
  
  pdf(paste0('/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/VolcanoPlots_', 
             variable_names[i], '_Cells_LCOnly.pdf'))
  print(tmp_graph)
  dev.off()
  
  png(paste0('/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/VolcanoPlots_', 
             variable_names[i], '_Cells_LCOnly.png'), width = 3000, height = 3000, res = 300)
  print(tmp_graph)
  dev.off()
  
  print(paste0('Plot done for ', variable_names[i], 
               if(use_break) ' (with axis break)' else ''))
}








































##### GSEA








# Load required libraries
library(fgsea)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(patchwork)  # Optional: for combining plots

# Define cell types
celltypes_vec <- c(#'All', 'PT', 'TAL', 'EC', 'DCTall', 'IC', 
  'POD', 'intercalated')

# Get gene sets from MSigDB for human
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

# Get GO gene sets
go_bp_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
go_bp_list <- split(go_bp_sets$gene_symbol, go_bp_sets$gs_name)

go_cc_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")
go_cc_list <- split(go_cc_sets$gene_symbol, go_cc_sets$gs_name)

go_mf_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")
go_mf_list <- split(go_mf_sets$gene_symbol, go_mf_sets$gs_name)

# Define gene set types
geneset_types <- list(
  Hallmark = hallmark_list,
  GO_BP = go_bp_list,
  GO_CC = go_cc_list,
  GO_MF = go_mf_list
)

# Create a list to store all results
all_gsea_results <- list()

# Loop through each cell type
for(celltype in celltypes_vec) {
  
  cat("\n=== Processing", celltype, "===\n")
  
  # Load your results file for this cell type
  de_results <- read.csv(paste0(dir.results, 'Full_NEBULA_', 
                                celltype, '_cells__LC_pooledoffset.csv'))
  
  # Create ranked gene list
  ranked_genes <- de_results %>%
    dplyr::select(logFC = summary.logFC_sexMale, 
                  gene_id = summary.gene) %>% 
    arrange(desc(logFC)) %>%
    pull(logFC, name = gene_id)
  
  # Loop through each gene set type
  for(geneset_name in names(geneset_types)) {
    
    cat("Running GSEA for", geneset_name, "...\n")
    
    # Run fGSEA
    gsea_results <- fgsea(
      pathways = geneset_types[[geneset_name]],
      stats = ranked_genes,
      minSize = 15,
      maxSize = 500,
      nperm = 10000
    )
    
    # Add cell type and gene set type columns
    gsea_results$celltype <- celltype
    gsea_results$geneset_type <- geneset_name
    
    # Store results
    result_name <- paste0(celltype, "_", geneset_name)
    all_gsea_results[[result_name]] <- gsea_results
    
    # Select top significant pathways
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
    
    # Only create plot if there are significant pathways
    if(nrow(top_pathways) > 0) {
      
      # Create dotplot
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
      
      # Save individual plot
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

# Combine all GSEA results into one dataframe
# Before saving, convert list columns to character strings
combined_gsea_save <- combined_gsea %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";")))

# Now save combined results
write.csv(combined_gsea_save, 
          paste0(dir.results, "GSEA/all_celltypes_all_genesets_gsea_results.csv"), 
          row.names = FALSE)

# Save separate files for each gene set type
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

# Summary of significant pathways per cell type and gene set
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














# Install required package if not already installed
if (!require("pdftools")) {
  install.packages("pdftools")
}

library(pdftools)

# Function to convert a single PDF to PNG
convert_pdf_to_png <- function(pdf_path, output_dir = NULL, dpi = 300) {
  # If no output directory specified, use same directory as PDF
  if (is.null(output_dir)) {
    output_dir <- dirname(pdf_path)
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get the base filename without extension
  base_name <- tools::file_path_sans_ext(basename(pdf_path))
  
  # Convert PDF to PNG
  # This creates one PNG per page
  png_files <- pdf_convert(
    pdf = pdf_path,
    format = "png",
    dpi = dpi,
    filenames = file.path(output_dir, paste0(base_name, "_page_%d.png"))
  )
  
  return(png_files)
}

# Convert all PDFs in a folder
convert_all_pdfs <- function(folder_path, output_dir = NULL, dpi = 300) {
  # Get all PDF files in the folder
  pdf_files <- list.files(folder_path, pattern = "\\.pdf$", 
                          full.names = TRUE, ignore.case = TRUE)
  
  if (length(pdf_files) == 0) {
    message("No PDF files found in the specified folder.")
    return(NULL)
  }
  
  message(paste("Found", length(pdf_files), "PDF file(s)"))
  
  # Convert each PDF
  all_png_files <- list()
  for (pdf in pdf_files) {
    message(paste("Converting:", basename(pdf)))
    png_files <- convert_pdf_to_png(pdf, output_dir, dpi)
    all_png_files[[basename(pdf)]] <- png_files
  }
  
  message("Conversion complete!")
  return(all_png_files)
}

# Example usage:
# Convert all PDFs in a folder
folder_path <- "/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/GSEA/"

# Or convert to a specific output directory
 converted_files <- convert_all_pdfs(folder_path, output_dir = "/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/GSEA/")

# Or convert with higher resolution
# converted_files <- convert_all_pdfs(folder_path, dpi = 600)









 
 #### Proteomics and Metabolomics 
 
 
 
 
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
 
 
 
 



 
 
 
 harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')
 
 dat <- harmonized_data %>% 
   dplyr::select(-dob) %>% 
   arrange(date_of_screen) %>% 
   dplyr::summarise(
     across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
     across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
     .by = c(record_id, visit)
   )
 
 dat <- dat %>% filter(group == 'Lean Control') %>% 
   filter(visit == 'baseline')
 
 
 data_dictionary <- readxl::read_xlsx('/Users/netio/Downloads/data_dictionary_master.xlsx')
 
 form_names <- unique(data_dictionary$form_name)
 proteo <- form_names[str_which(form_names, pattern = 'proteom')]
 metab <- form_names[str_which(form_names, pattern = 'metab')]
 
 variables_class <- c('proteomics', 'metabolomics', 'metabolomics_blood_raw',
                      'az_urine_metabolites', 'metabolomics_aq')
 
 data_dictionary_small <- data_dictionary %>% 
   filter(form_name %in% variables_class)
 
 
 
 
 
 # Fix data types before creating the table
 library(gtsummary)
 library(gt)
 library(dplyr)
 
 # Convert variables to proper data types
 combined_df <- dat %>% 
   filter(!is.na(citrate)) %>% 
   mutate(
     # Ensure continuous variables are numeric
     age = as.numeric(age),
     bmi = as.numeric(bmi),
     hba1c = as.numeric(hba1c),
     
     # Ensure categorical variables are factors or characters
     sex = as.factor(sex),
     race_ethnicity = as.factor(race_ethnicity),
     study = as.factor(study),
     group = as.factor(group)
   )
 
 
 
 # Now create the table with proper data types
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
     # Skip categorical p-values if they cause issues
   )) %>%
   add_overall(col_label = "**Overall**\nN = {N}") %>%
   modify_header(label ~ "**Characteristic**") %>%
   modify_spanning_header(all_stat_cols() ~ "**Sex**") %>%
   modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")
 
 # Save version with epic
 desc_table1_fixed %>%
   as_gt() %>%
   tab_options(
     table.font.size = 11,
     heading.title.font.size = 14,
     column_labels.font.size = 12
   ) %>%
   gtsave("C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/omics/LeanControl_Omics_Demographics.png", 
          vwidth = 1200, vheight = 800)
 
 
 
 
 
 
 
 
 
 
 
 
 
 # Find which columns actually exist
 existing_cols <- intersect(data_dictionary_small$variable_name, names(dat))
 
 # Select only those
 dat_omics <- dat %>% 
   dplyr::select(record_id, group, sex, study, all_of(existing_cols))
 
 # Check how many were found vs missing
 cat("Found:", length(existing_cols), "out of", length(data_dictionary_small$variable_name), "\n")
 cat("Missing:", length(data_dictionary_small$variable_name) - length(existing_cols), "\n")
 
 
 
 ##Analysis
 
 library(dplyr)
 library(ggplot2)
 library(tidyr)
 library(gridExtra)
 library(readxl)
 
 # Read data dictionary
 data_dictionary <- readxl::read_xlsx('/Users/netio/Downloads/data_dictionary_master.xlsx')
 
 # Set output folder
 output_folder <- '/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/omics/'
 
 # Get all numeric columns (excluding ID columns and grouping variables)
 numeric_cols <- names(dat_omics)[!names(dat_omics) %in% c('record_id', 'group', 'sex')]
 
 # Initialize results dataframe
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
 
 # Loop through each variable
 for(var in numeric_cols){
   
   # Get label and form_name from data dictionary
   var_label <- data_dictionary$label[data_dictionary$variable_name == var]
   var_form <- data_dictionary$form_name[data_dictionary$variable_name == var]
   
   if(length(var_label) == 0) var_label <- var  # Use variable name if no label found
   if(length(var_form) == 0) var_form <- NA  # Use NA if no form_name found
   
   # Create temporary dataframe
   temp_data <- dat_omics %>%
     filter(!is.na(sex) & sex != "") %>%
     select(sex, all_of(var)) %>%
     rename(value = all_of(var)) %>%
     filter(!is.na(value))
   
   # Skip if not enough data
   if(nrow(temp_data) < 6 || length(unique(temp_data$sex)) < 2) next
   
   # Split by sex
   male_data <- temp_data$value[temp_data$sex == "Male"]
   female_data <- temp_data$value[temp_data$sex == "Female"]
   
   # Skip if either group has less than 3 observations
   if(length(male_data) < 3 || length(female_data) < 3) next
   
   # Test for normality (if sample size allows)
   normal_male <- if(length(male_data) >= 3 & length(male_data) <= 5000) {
     shapiro.test(male_data)$p.value > 0.05
   } else TRUE
   
   normal_female <- if(length(female_data) >= 3 & length(female_data) <= 5000) {
     shapiro.test(female_data)$p.value > 0.05
   } else TRUE
   
   # Choose appropriate test
   if(normal_male & normal_female){
     # Use t-test
     test_result <- t.test(value ~ sex, data = temp_data, var.equal = FALSE)
     test_name <- "Welch's t-test"
   } else {
     # Use Mann-Whitney U test
     test_result <- wilcox.test(value ~ sex, data = temp_data)
     test_name <- "Mann-Whitney U"
   }
   
   # Store results
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
 
 # Add corrected p-values (multiple methods)
 results <- results %>%
   mutate(
     p_bonferroni = p.adjust(p_value, method = "bonferroni"),
     p_fdr = p.adjust(p_value, method = "fdr"),  # Benjamini-Hochberg
     difference = mean_male - mean_female,
     abs_difference = abs(difference)
   ) %>%
   arrange(p_value)
 
 # View top results
 print(head(results, 20))
 
 # Save full results
 write.csv(results, paste0(output_folder, 'sex_comparison_results.csv'), row.names = FALSE)
 
 # Plot top 10 most significant
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
   
   # Get p-value for this variable
   p_val <- top_results$p_value[i]
   p_fdr <- top_results$p_fdr[i]
   test_type <- top_results$test_used[i]
   
   # Format p-values
   p_label <- ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val))
   p_fdr_label <- ifelse(p_fdr < 0.001, "FDR < 0.001", sprintf("FDR = %.3f", p_fdr))
   
   plot_list[[i]] <- ggplot(plot_data, aes(x = sex, y = value, fill = sex)) +
     geom_boxplot() +
     scale_fill_manual(values = c("Male" = "#00008B", "Female" = "#8B0000")) +
     theme_classic(base_size = 12) +
     labs(
       title = var_label,  # Use label instead of variable code
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
 
 # Save plots
 png(paste0(output_folder, 'top10_sex_differences.png'), 
     height = 15, width = 12, units = 'in', res = 300)
 grid.arrange(grobs = plot_list, ncol = 2)
 dev.off()
 
 # Also create a plot for top 20
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
       title = var_label,  # Use label instead of variable code
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
 
 # Summary by form_name
 print("\nSummary by analysis method (form_name):")
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
 
 # Prepare your data
 # Assume you have a dataframe with metabolite names and statistics
 # metabolites_data should have columns: metabolite_name, fold_change, p_value
 
 # Option A: Over-representation analysis (ORA)
 # Use if you have a list of significantly changed metabolites
 
 # Create a vector of significant metabolites
 
 metabolites_data <- data.table::fread('/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/omics/sex_comparison_results.csv')
   
   
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
 mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway", 2) # or "kegg_pathway"
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
 # Use if you have quantitative data (fold changes, concentrations)
 
 # Prepare data with fold changes
 metabolite_scores <- metabolites_data %>%
   select(metabolite_name, fold_change) %>%
   na.omit()
 
 # Write to file (MetaboAnalyst format)
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
 
 
 # ===== METHOD 3: Using FGSEA (Gene Set Enrichment Analysis adapted for metabolomics) =====
 library(fgsea)
 
 # You'll need pathway databases
 # Example: Load KEGG pathways (you need to create/download this)
 # Format: list where each element is a pathway with metabolite IDs
 
 # Load or create pathway database
 # kegg_pathways <- list(
 #   "Glycolysis" = c("HMDB0000122", "HMDB0000094", ...),
 #   "TCA Cycle" = c("HMDB0000094", "HMDB0000156", ...),
 #   ...
 # )
 
 # Prepare ranked list of metabolites
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
 
 # View results
 fgsea_results <- fgsea_results %>%
   arrange(padj)
 
 print(head(fgsea_results, 20))
 
 # Plot enrichment
 plotEnrichment(kegg_pathways[["Glycolysis"]], metabolite_ranks) +
   labs(title = "Glycolysis Pathway Enrichment")
 
 
 # ===== METHOD 4: Manual pathway enrichment using hypergeometric test =====
 
 # Load pathway database (example structure)
 pathway_db <- list(
   "Glycolysis / Gluconeogenesis" = c("Glucose", "Pyruvate", "Lactate", "Glucose-6-phosphate"),
   "TCA Cycle" = c("Citrate", "Succinate", "Fumarate", "Malate"),
   "Fatty Acid Metabolism" = c("Palmitate", "Stearate", "Oleate", "Acetyl-CoA")
   # ... add more pathways
 )
 
 # Your significant metabolites
 sig_metabolites <- metabolites_data %>%
   filter(p_value < 0.05) %>%
   pull(metabolite_name)
 
 # Total metabolites tested
 total_metabolites <- nrow(metabolites_data)
 
 # Perform enrichment for each pathway
 enrichment_results <- data.frame()
 
 for(pathway_name in names(pathway_db)) {
   pathway_metabolites <- pathway_db[[pathway_name]]
   
   # Overlap between significant and pathway
   overlap <- length(intersect(sig_metabolites, pathway_metabolites))
   
   # Hypergeometric test
   p_val <- phyper(
     q = overlap - 1,
     m = length(pathway_metabolites),  # metabolites in pathway
     n = total_metabolites - length(pathway_metabolites),  # metabolites not in pathway
     k = length(sig_metabolites),  # significant metabolites
     lower.tail = FALSE
   )
   
   enrichment_results <- rbind(enrichment_results, data.frame(
     pathway = pathway_name,
     overlap = overlap,
     pathway_size = length(pathway_metabolites),
     p_value = p_val
   ))
 }
 
 # Adjust for multiple testing
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
 
 
 # ===== METHOD 5: Using online MetaboAnalyst (recommended for beginners) =====
 # Prepare data for upload to https://www.metaboanalyst.ca/
 
 # Format 1: For ORA (just list of metabolites)
 write.table(significant_metabolites,
             file = paste0(dir.results, "significant_metabolites_for_metaboanalyst.txt"),
             row.names = FALSE, col.names = FALSE, quote = FALSE)
 
 # Format 2: For MSEA (metabolites with scores)
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
   
   # Remove NAs and duplicates
   unique_ids <- unique(uniprot_ids[!is.na(uniprot_ids)])
   
   if(length(unique_ids) == 0) {
     return(data.frame(Entry = character(), 
                       Gene.Names = character(), 
                       Protein.names = character()))
   }
   
   # Split into batches of 100 (API limit)
   batch_size <- 100
   batches <- split(unique_ids, ceiling(seq_along(unique_ids) / batch_size))
   
   all_results <- data.frame()
   
   cat("Retrieving protein names from UniProt...\n")
   
   for (i in seq_along(batches)) {
     cat(paste0("Processing batch ", i, " of ", length(batches), "...\n"))
     
     batch <- batches[[i]]
     
     # Query UniProt API
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
       # Parse response
       content <- content(response, "text", encoding = "UTF-8")
       if(nchar(content) > 0) {
         batch_results <- read.delim(text = content, sep = "\t", stringsAsFactors = FALSE)
         all_results <- rbind(all_results, batch_results)
       }
     } else {
       warning(paste("Failed to retrieve batch", i, ". Status code:", status_code(response)))
     }
     
     # Be nice to the API
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
 
 dat <- dat %>% filter(group == 'Lean Control')
 
 # Read data dictionary
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
 
 # Extract UniProt IDs (pattern matches standard UniProt accession format)
 # This regex matches patterns like P12345, Q9Y123, O95238, etc.
 uniprot_ids <- proteomics_vars$label %>%
   str_extract("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")
 
 # Get protein names from UniProt
 protein_mapping <- get_protein_names(uniprot_ids)
 
 # Create enhanced data dictionary with protein information
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
     # Create clean gene name (take first gene if multiple)
     gene_name_clean = ifelse(!is.na(gene_name), 
                              str_trim(str_split(gene_name, " ", simplify = TRUE)[,1]),
                              NA),
     # Create enhanced label with protein info
     label_enhanced = case_when(
       !is.na(gene_name_clean) & !is.na(protein_name) ~ paste0(gene_name_clean, " - ", protein_name),
       !is.na(gene_name_clean) ~ gene_name_clean,
       TRUE ~ label
     )
   )
 
 cat(paste0("\nSuccessfully mapped ", sum(!is.na(data_dictionary_enhanced$gene_name)), 
            " proteins out of ", nrow(proteomics_vars), " proteomics variables\n\n"))
 
 # Replace original data dictionary with enhanced version
 data_dictionary <- data_dictionary_enhanced
 
 data_dictionary_small <- data_dictionary %>% 
   filter(form_name %in% variables_class)
 
 data_dictionary_small$variable_name <- str_replace_all(data_dictionary_small$variable_name, pattern = '_', replacement = '.')
 
 # Find which columns actually exist
 existing_cols <- intersect(data_dictionary_small$variable_name, names(dat))
 
 # Select only those
 dat_omics <- dat %>% 
   dplyr::select(record_id, group, study, sex, age, all_of(existing_cols))
 
 # Check how many were found vs missing
 cat("Found:", length(existing_cols), "out of", length(data_dictionary_small$variable_name), "\n")
 cat("Missing:", length(data_dictionary_small$variable_name) - length(existing_cols), "\n")
 
 ##Analysis
 
 # Set output folder
 output_folder <- '/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/omics/'
 
 # Get all numeric columns (excluding ID columns and grouping variables)
 numeric_cols <- names(dat_omics)[!names(dat_omics) %in% c('record_id', 'group', 'sex')]
 
 # Initialize results dataframe
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
 
 # Loop through each variable
 for(var in numeric_cols){
   
   # Get all information from enhanced data dictionary
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
   
   # Create temporary dataframe
   temp_data <- dat_omics %>%
     filter(!is.na(sex) & sex != "") %>%
     select(sex, all_of(var)) %>%
     rename(value = all_of(var)) %>%
     filter(!is.na(value))
   
   # Skip if not enough data
   if(nrow(temp_data) < 6 || length(unique(temp_data$sex)) < 2) next
   
   # Split by sex
   male_data <- temp_data$value[temp_data$sex == "Male"]
   female_data <- temp_data$value[temp_data$sex == "Female"]
   
   # Skip if either group has less than 3 observations
   if(length(male_data) < 3 || length(female_data) < 3) next
   
   # Test for normality (if sample size allows)
   normal_male <- if(length(male_data) >= 3 & length(male_data) <= 5000) {
     shapiro.test(male_data)$p.value > 0.05
   } else TRUE
   
   normal_female <- if(length(female_data) >= 3 & length(female_data) <= 5000) {
     shapiro.test(female_data)$p.value > 0.05
   } else TRUE
   
   # Choose appropriate test
   if(normal_male & normal_female){
     # Use t-test
     test_result <- t.test(value ~ sex, data = temp_data, var.equal = FALSE)
     test_name <- "Welch's t-test"
   } else {
     # Use Mann-Whitney U test
     test_result <- wilcox.test(value ~ sex, data = temp_data)
     test_name <- "Mann-Whitney U"
   }
   
   # Store results
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
 
 # Add corrected p-values (multiple methods)
 results <- results %>%
   mutate(
     p_bonferroni = p.adjust(p_value, method = "bonferroni"),
     p_fdr = p.adjust(p_value, method = "fdr"),  # Benjamini-Hochberg
     difference = mean_male - mean_female,
     abs_difference = abs(difference)
   ) %>%
   arrange(p_value)
 
 # View top results
 cat("\n=== Top 20 Results ===\n")
 print(head(results %>% select(gene_name, protein_name, p_value, p_fdr, mean_male, mean_female), 20))
 
 # Save full results
 write.csv(results, paste0(output_folder, 'proteomics_sex_comparison_results.csv'), row.names = FALSE)
 
 cat("\n=== Creating plots ===\n")
 
 # Plot top 10 most significant
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
   
   # Get p-value for this variable
   p_val <- top_results$p_value[i]
   p_fdr <- top_results$p_fdr[i]
   test_type <- top_results$test_used[i]
   
   # Format p-values
   p_label <- ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val))
   p_fdr_label <- ifelse(p_fdr < 0.001, "FDR < 0.001", sprintf("FDR = %.3f", p_fdr))
   
   plot_list[[i]] <- ggplot(plot_data, aes(x = sex, y = value, fill = sex)) +
     geom_boxplot() +
     scale_fill_manual(values = c("Male" = "#00008B", "Female" = "#8B0000")) +
     theme_classic(base_size = 12) +
     labs(
       title = var_label,  # Now uses gene name and protein name
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
 
 # Save plots
 png(paste0(output_folder, 'proteomics_top10_sex_differences.png'), 
     height = 15, width = 12, units = 'in', res = 300)
 grid.arrange(grobs = plot_list, ncol = 2)
 dev.off()
 
 # Also create a plot for top 20
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
       title = var_label,  # Now uses gene name and protein name
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
 
 # Summary by form_name
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
 
 # Summary of proteins with gene names mapped
 cat("\nProteomics mapping summary:\n")
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
 
 # Function to get protein names from UniProt IDs
 get_protein_names <- function(uniprot_ids) {
   
   # Remove NAs and duplicates
   unique_ids <- unique(uniprot_ids[!is.na(uniprot_ids)])
   
   if(length(unique_ids) == 0) {
     return(data.frame(Entry = character(), 
                       Gene.Names = character(), 
                       Protein.names = character()))
   }
   
   # Split into batches of 100 (API limit)
   batch_size <- 100
   batches <- split(unique_ids, ceiling(seq_along(unique_ids) / batch_size))
   
   all_results <- data.frame()
   
   cat("Retrieving protein names from UniProt...\n")
   
   for (i in seq_along(batches)) {
     cat(paste0("Processing batch ", i, " of ", length(batches), "...\n"))
     
     batch <- batches[[i]]
     
     # Query UniProt API
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
       # Parse response
       content <- content(response, "text", encoding = "UTF-8")
       if(nchar(content) > 0) {
         batch_results <- read.delim(text = content, sep = "\t", stringsAsFactors = FALSE)
         all_results <- rbind(all_results, batch_results)
       }
     } else {
       warning(paste("Failed to retrieve batch", i, ". Status code:", status_code(response)))
     }
     
     # Be nice to the API
     Sys.sleep(0.5)
   }
   
   cat(paste0("Retrieved information for ", nrow(all_results), " proteins\n"))
   
   return(all_results)
 }
 
 # Read your existing results
 results <- read.csv('/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/omics/proteomics_sex_comparison_results.csv')
 
 cat("Original results dimensions:", nrow(results), "rows\n\n")
 
 # Read data dictionary
 data_dictionary <- readxl::read_xlsx('/Users/netio/Downloads/data_dictionary_master.xlsx')
 
 # Filter to proteomics only
 proteomics_dict <- data_dictionary %>%
   filter(form_name == 'proteomics')
 
 cat("Proteomics variables in dictionary:", nrow(proteomics_dict), "\n")
 cat("\nFirst few labels from dictionary:\n")
 print(head(proteomics_dict$label, 10))
 
 # Extract UniProt IDs from the label column
 # Pattern matches: P12345, Q9Y123, O95238, etc.
 proteomics_dict <- proteomics_dict %>%
   mutate(
     uniprot_id = str_extract(label, "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")
   )
 
 # Check how many UniProt IDs we found
 cat("\nUniProt IDs found:", sum(!is.na(proteomics_dict$uniprot_id)), "\n")
 cat("Sample UniProt IDs:\n")
 print(head(proteomics_dict$uniprot_id[!is.na(proteomics_dict$uniprot_id)], 10))
 
 # Get unique UniProt IDs
 uniprot_ids <- unique(proteomics_dict$uniprot_id)
 uniprot_ids <- uniprot_ids[!is.na(uniprot_ids)]
 
 cat("\nUnique UniProt IDs to query:", length(uniprot_ids), "\n\n")
 
 # Get protein names from UniProt
 protein_mapping <- get_protein_names(uniprot_ids)
 
 # Clean up the protein mapping
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
 
 # Create enhanced dictionary with protein info
 proteomics_dict_enhanced <- proteomics_dict %>%
   left_join(protein_mapping, by = "uniprot_id") %>%
   mutate(
     # Create display label: Gene Name - Protein Description
     label_enhanced = case_when(
       !is.na(gene_name_clean) & !is.na(protein_name) ~ paste0(gene_name_clean, " - ", protein_name),
       !is.na(gene_name_clean) ~ gene_name_clean,
       TRUE ~ label
     )
   )
 
 # Now merge with results
 # First, need to match variable names (results has dots, dictionary might have underscores)
 results_enhanced <- results %>%
   mutate(variable_underscore = str_replace_all(variable, "\\.", "_")) %>%
   left_join(
     proteomics_dict_enhanced %>% 
       select(variable_name, uniprot_id, gene_name_clean, gene_name_full, protein_name, label_enhanced),
     by = c("variable_underscore" = "variable_name")
   ) %>%
   mutate(
     # Update columns with new information
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
 
 # Check mapping success
 cat("\n=== Mapping Summary ===\n")
 cat("Total results:", nrow(results_enhanced), "\n")
 cat("Results with gene names:", sum(!is.na(results_enhanced$gene_name) & results_enhanced$gene_name != ""), "\n")
 cat("Results with protein names:", sum(!is.na(results_enhanced$protein_name) & results_enhanced$protein_name != ""), "\n")
 
 # View top results
 cat("\n=== Top 10 Results with Protein Names ===\n")
 print(results_enhanced %>% 
         filter(!is.na(gene_name)) %>%
         select(gene_name, protein_name, p_value, p_fdr, mean_male, mean_female) %>%
         head(10))
 
 # Save enhanced results
 output_file <- '/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/omics/proteomics_sex_comparison_results_ENHANCED.csv'
 write.csv(results_enhanced, output_file, row.names = FALSE)
 
 cat("\n=== Saved enhanced results to: ===\n")
 cat(output_file, "\n")
 
 cat("\n=== Complete! ===\n")
 cat("New columns added:\n")
 cat("  - variable_label_enhanced: Gene name and protein description\n")
 cat("  - uniprot_id: UniProt accession number\n")
 cat("  - gene_name: Primary gene name\n")
 cat("  - gene_name_full: All gene names/aliases\n")
 cat("  - protein_name: Full protein description\n") 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ########## Lineage Tracing (Slingshot) Analysis in Lean Controls (PT Cells)
 
 library(slingshot)
 #library(condiments)
 
 
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
 
 
 
 
 
 
 
 load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')
 
 so_subset <- so_kpmp_sc
 remove(so_kpmp_sc)
 
 #dat_groups <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_GroupAssignments.txt')
 #dat_groups <- dat_groups %>% filter(group2 %in% c('Lean Control', 'T2D-No SGLTi2'))
 
 #so_subset <- subset(so_subset, record_id == dat_groups$record_id)
 test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))
 
 
 
 
 
 
 
 
 dir.results <- 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/pseudotime/'
 
 
 so_subset <- subset(so_subset, subset = record_id != 'CRC-55')
 so_subset <- subset(so_subset, subset = group == 'Lean_Control')
 test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))
 
 
 
 #PT Cells
 so_subset <- subset(so_subset, subset = celltype2 == 'PT')
 so_subset <- RunUMAP(so_subset, dims = 1:30)
 
 sling_res <- slingshot(as.SingleCellExperiment(so_subset), clusterLabels = 'KPMP_celltype', 
                        start.clus = 'PT-S1', end.clus = 'aPT', reducedDim = 'UMAP')
 
 so_subset$pseudotime <- slingPseudotime(sling_res)[,1]
 
 
 
 
 
 
 library(ggplot2)
 library(viridis)
 
 # 1. Plot slingshot lineage tracing on UMAP
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
   sex = so_subset$sex  # Assuming you have a 'sex' column in metadata
 )
 
 # Plot 1: UMAP with pseudotime
 p1 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
   geom_point(size = 0.5, alpha = 0.6) +
   scale_color_viridis(option = "plasma", na.value = "grey80") +
   theme_classic() +
   labs(title = "Slingshot Pseudotime on UMAP",
        color = "Pseudotime") +
   theme(legend.position = "right")
 
 # Add slingshot curves
 for(i in seq_along(sling_curves)) {
   curve_coords <- sling_curves[[i]]$s[sling_curves[[i]]$ord, ]
   p1 <- p1 + geom_path(data = data.frame(UMAP_1 = curve_coords[, 1], 
                                          UMAP_2 = curve_coords[, 2]),
                        aes(x = UMAP_1, y = UMAP_2),
                        color = "black", size = 2, inherit.aes = FALSE)
 }
 
 print(p1)
 ggsave(paste0(dir.results, "Slingshot_UMAP_Pseudotime.pdf"), 
        plot = p1, width = 8, height = 6)
 
 # Plot 2: UMAP with cell types
 p2 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
   geom_point(size = 0.5, alpha = 0.6) +
   theme_classic() +
   labs(title = "Cell Types with Slingshot Trajectory",
        color = "Cell Type") +
   theme(legend.position = "right")
 
 # Add slingshot curves
 for(i in seq_along(sling_curves)) {
   curve_coords <- sling_curves[[i]]$s[sling_curves[[i]]$ord, ]
   p2 <- p2 + geom_path(data = data.frame(UMAP_1 = curve_coords[, 1], 
                                          UMAP_2 = curve_coords[, 2]),
                        aes(x = UMAP_1, y = UMAP_2),
                        color = "black", size = 2, inherit.aes = FALSE)
 }
 
 print(p2)
 ggsave(paste0(dir.results, "Slingshot_UMAP_Celltypes.pdf"), 
        plot = p2, width = 8, height = 6)
 
 # 2. Compare pseudotime between males and females
 # Remove cells with NA pseudotime
 plot_df_clean <- plot_df %>% filter(!is.na(pseudotime))
 
 # Summary statistics
 pseudotime_summary <- plot_df_clean %>%
   group_by(sex) %>%
   summarise(
     n = n(),
     mean_pseudotime = mean(pseudotime),
     median_pseudotime = median(pseudotime),
     sd_pseudotime = sd(pseudotime),
     se_pseudotime = sd(pseudotime) / sqrt(n())
   )
 
 print(pseudotime_summary)
 write.csv(pseudotime_summary, paste0(dir.results, "Pseudotime_Summary_by_Sex.csv"), 
           row.names = FALSE)
 
 # Statistical test - Wilcoxon test
 stat_test <- plot_df_clean %>%
   wilcox_test(pseudotime ~ sex) %>%
   add_significance()
 
 print(stat_test)
 
 # Plot 3: Violin plot comparing pseudotime by sex
 p3 <- ggplot(plot_df_clean, aes(x = sex, y = pseudotime, fill = sex)) +
   geom_violin(alpha = 0.6, trim = FALSE) +
   geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
   geom_jitter(width = 0.1, alpha = 0.1, size = 0.5) +
   stat_compare_means(method = "wilcox.test", label = "p.format") +
   theme_classic() +
   labs(title = "Pseudotime Comparison Between Sexes",
        x = "Sex",
        y = "Pseudotime") +
   theme(legend.position = "none")
 
 print(p3)
 ggsave(paste0(dir.results, "Pseudotime_by_Sex_Violin.pdf"), 
        plot = p3, width = 6, height = 6)
 
 # Plot 4: Density plot of pseudotime distribution by sex
 p4 <- ggplot(plot_df_clean, aes(x = pseudotime, fill = sex)) +
   geom_density(alpha = 0.5) +
   theme_classic() +
   labs(title = "Pseudotime Distribution by Sex",
        x = "Pseudotime",
        y = "Density",
        fill = "Sex")
 
 print(p4)
 ggsave(paste0(dir.results, "Pseudotime_by_Sex_Density.pdf"), 
        plot = p4, width = 8, height = 6)
 
 # NEW: Plot 4b: Density plot over pseudotime trajectory (all cells)
 p4b <- ggplot(plot_df_clean, aes(x = pseudotime)) +
   geom_density(fill = "steelblue", alpha = 0.6, color = "black") +
   theme_classic() +
   labs(title = "Cell Density Along Pseudotime Trajectory",
        x = "Pseudotime",
        y = "Density") +
   theme(plot.title = element_text(hjust = 0.5))
 
 print(p4b)
 ggsave(paste0(dir.results, "Pseudotime_Density_Overall.pdf"), 
        plot = p4b, width = 8, height = 5)
 
 # NEW: Plot 4c: Stacked density plot showing progression by sex
 p4c <- ggplot(plot_df_clean, aes(x = pseudotime, fill = sex, color = sex)) +
   geom_density(alpha = 0.3, size = 1) +
   theme_classic() +
   labs(title = "Cell Density Along Pseudotime Trajectory by Sex",
        x = "Pseudotime",
        y = "Density",
        fill = "Sex",
        color = "Sex") +
   theme(plot.title = element_text(hjust = 0.5))
 
 print(p4c)
 ggsave(paste0(dir.results, "Pseudotime_Density_by_Sex.pdf"), 
        plot = p4c, width = 8, height = 5)
 
 # NEW: Plot 4d: Ridgeline-style density plot by cell type along pseudotime
 library(ggridges)
 
 p4d <- ggplot(plot_df_clean, aes(x = pseudotime, y = celltype, fill = celltype)) +
   geom_density_ridges(alpha = 0.6, scale = 2) +
   theme_classic() +
   labs(title = "Cell Type Distribution Along Pseudotime",
        x = "Pseudotime",
        y = "Cell Type") +
   theme(legend.position = "none",
         plot.title = element_text(hjust = 0.5))
 
 print(p4d)
 ggsave(paste0(dir.results, "Pseudotime_Density_by_Celltype_Ridge.pdf"), 
        plot = p4d, width = 8, height = 6)
 
 # Plot 5: UMAP split by sex
 p5 <- ggplot(plot_df_clean, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
   geom_point(size = 0.5, alpha = 0.6) +
   scale_color_viridis(option = "plasma") +
   facet_wrap(~sex) +
   theme_classic() +
   labs(title = "Slingshot Pseudotime by Sex",
        color = "Pseudotime")
 
 # Add curves to each facet
 for(i in seq_along(sling_curves)) {
   curve_coords <- sling_curves[[i]]$s[sling_curves[[i]]$ord, ]
   p5 <- p5 + geom_path(data = data.frame(UMAP_1 = curve_coords[, 1], 
                                          UMAP_2 = curve_coords[, 2]),
                        aes(x = UMAP_1, y = UMAP_2),
                        color = "black", size = 1.5, inherit.aes = FALSE)
 }
 
 print(p5)
 ggsave(paste0(dir.results, "Slingshot_UMAP_by_Sex.pdf"), 
        plot = p5, width = 12, height = 5)
 
 # Combined panel plot - updated to include density over pseudotime
 combined_plot <- (p1 | p2) / (p4c | p3)
 print(combined_plot)
 ggsave(paste0(dir.results, "Slingshot_Analysis_Combined.pdf"), 
        plot = combined_plot, width = 14, height = 10)
 
 # Alternative combined plot with ridgeline
 combined_plot2 <- (p1 | p2) / (p4c | p4d)
 print(combined_plot2)
 ggsave(paste0(dir.results, "Slingshot_Analysis_Combined_with_Ridge.pdf"), 
        plot = combined_plot2, width = 14, height = 10)
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
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
   celltype = so_subset$KPMP_celltype,  # Using original cell type labels
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
 
 # Define the three regions
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
   # Extract cells in this region
   region_cells <- plot_df_clean %>%
     filter(pseudotime >= region$start & pseudotime <= region$end)
   
   # Calculate counts in this region BY SEX
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
   
   # Add any missing cell type/sex combinations with 0%
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
 
 # Print detailed summary
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
 
 # Save the detailed composition data
 write.csv(all_region_data, 
           paste0(dir.results, "Pseudotime_Regions_Celltype_Sex_Percentage.csv"), 
           row.names = FALSE)
 
 # Create bar plots showing percentage of each cell type in each region, SPLIT BY SEX
 celltype_colors <- c("PT-S1/S2" = "#4DAF4A",  # Green
                      "PT-S3" = "#377EB8",      # Blue
                      "aPT" = "#E41A1C")        # Red
 
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
     ylim(0, max(all_region_data$percent_of_celltype) * 1.15)  # Add space for labels
 }
 
 # Plot G - Ridge plot showing pseudotime distribution by cell type
 # Order cell types by median pseudotime
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
   # Add vertical lines for region boundaries
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
 
 # Label bar charts
 for(i in 1:length(bar_charts)) {
   bar_charts[[i]] <- bar_charts[[i]] + 
     labs(tag = LETTERS[3 + i]) +
     theme(plot.tag = element_text(size = 16, face = "bold"))
 }
 
 # Label ridge plot as G
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
 
 
 # Create a summary table showing key differences
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
 #library(condiments)
 
 
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
 
 
 
 
 
 load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')
 
 so_subset <- so_kpmp_sc
 remove(so_kpmp_sc)
 
 #dat_groups <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_GroupAssignments.txt')
 #dat_groups <- dat_groups %>% filter(group2 %in% c('Lean Control', 'T2D-No SGLTi2'))
 
 #so_subset <- subset(so_subset, record_id == dat_groups$record_id)
 test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))
 
 
 
 
 
 
 
 
 dir.results <- 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/pseudotime/'
 
 
 so_subset <- subset(so_subset, subset = record_id != 'CRC-55')
 so_subset <- subset(so_subset, subset = group == 'Lean_Control')
 test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))
 
 
 
 #PT Cells
 so_subset <- subset(so_subset, subset = celltype2 == 'PT')
 so_subset <- RunUMAP(so_subset, dims = 1:30)
 
 sling_res <- slingshot(as.SingleCellExperiment(so_subset), clusterLabels = 'KPMP_celltype', 
                        start.clus = 'PT-S1', end.clus = 'aPT', reducedDim = 'UMAP')
 
 so_subset$pseudotime <- slingPseudotime(sling_res)[,1]
 
 
 
 
 
 
 library(ggplot2)
 library(viridis)
 
 # 1. Plot slingshot lineage tracing on UMAP
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
   sex = so_subset$sex  # Assuming you have a 'sex' column in metadata
 )
 
 
 
 setwd('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/LeanControl_Only/pseudotime/further_exploration/')
 
 
 library(tradeSeq)
 library(phenopath)
 library(dplyr)
 library(ggplot2)
 library(viridis)
 library(Seurat)
 library(patchwork)
 library(ggridges)
 library(effsize)
 library(emmeans)
 library(clusterProfiler)
 library(org.Hs.eg.db)
 library(VennDiagram)
 library(grid)
 
 # ============================================
 # DATA VALIDATION AND CLEANING
 # ============================================
 
 print("Starting comprehensive sex-based trajectory analysis...")
 print("========================================")
 
 # Check for missing values in key variables
 print("Checking for missing values:")
 print(paste("Missing pseudotime:", sum(is.na(plot_df$pseudotime))))
 print(paste("Missing sex:", sum(is.na(plot_df$sex))))
 print(paste("Missing UMAP_1:", sum(is.na(plot_df$UMAP_1))))
 print(paste("Missing UMAP_2:", sum(is.na(plot_df$UMAP_2))))
 
 # Remove rows with missing pseudotime or sex
 plot_df_clean <- plot_df %>%
   filter(!is.na(pseudotime) & !is.na(sex) & !is.na(UMAP_1) & !is.na(UMAP_2))
 
 print(paste("Original rows:", nrow(plot_df)))
 print(paste("Clean rows:", nrow(plot_df_clean)))
 print(paste("Removed rows:", nrow(plot_df) - nrow(plot_df_clean)))
 
 # Check sex distribution
 print("Sex distribution:")
 print(table(plot_df_clean$sex))
 
 # Make sure sex is a factor with exactly 2 levels
 plot_df_clean$sex <- factor(plot_df_clean$sex)
 if(length(levels(plot_df_clean$sex)) != 2) {
   stop("Sex variable must have exactly 2 levels")
 }
 
 # ============================================
 # PART 1: Test for Overall Trajectory Differences
 # ============================================
 
 print("\n========================================")
 print("PART 1: OVERALL TRAJECTORY ANALYSIS")
 print("========================================\n")
 
 # 1. Compare pseudotime distributions between sexes
 trajectory_stats <- plot_df_clean %>%
   group_by(sex) %>%
   summarise(
     n_cells = n(),
     mean_pt = mean(pseudotime, na.rm = TRUE),
     median_pt = median(pseudotime, na.rm = TRUE),
     sd_pt = sd(pseudotime, na.rm = TRUE),
     min_pt = min(pseudotime, na.rm = TRUE),
     max_pt = max(pseudotime, na.rm = TRUE),
     q25_pt = quantile(pseudotime, 0.25, na.rm = TRUE),
     q75_pt = quantile(pseudotime, 0.75, na.rm = TRUE)
   )
 
 print("Trajectory Statistics by Sex:")
 print(trajectory_stats)
 
 # 2. Statistical tests for pseudotime differences
 # Wilcoxon rank-sum test
 wilcox_result <- wilcox.test(pseudotime ~ sex, data = plot_df_clean)
 print(paste("\nWilcoxon test p-value:", wilcox_result$p.value))
 
 # T-test (if distributions are relatively normal)
 t_result <- t.test(pseudotime ~ sex, data = plot_df_clean)
 print(paste("T-test p-value:", t_result$p.value))
 
 # Effect size (Cohen's d)
 cohens_d <- cohen.d(plot_df_clean$pseudotime, plot_df_clean$sex)
 print(paste("Cohen's d effect size:", cohens_d$estimate))
 
 # 3. Visualize trajectory differences
 p1 <- ggplot(plot_df_clean, aes(x = pseudotime, y = sex, fill = sex)) +
   geom_density_ridges(alpha = 0.7, scale = 0.9) +
   labs(title = "Pseudotime Distribution by Sex",
        x = "Pseudotime",
        y = "Sex") +
   theme_ridges() +
   theme(legend.position = "none")
 
 p2 <- ggplot(plot_df_clean, aes(x = sex, y = pseudotime, fill = sex)) +
   geom_violin(alpha = 0.7) +
   geom_boxplot(width = 0.2, alpha = 0.5) +
   labs(title = "Pseudotime by Sex",
        x = "Sex",
        y = "Pseudotime") +
   theme_classic() +
   theme(legend.position = "none")
 
 p3 <- ggplot(plot_df_clean, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
   geom_point(size = 0.5, alpha = 0.6) +
   facet_wrap(~sex) +
   scale_color_viridis(option = "magma") +
   labs(title = "Trajectory by Sex on UMAP") +
   theme_classic()
 
 # Combine plots
 print((p1 / p2) | p3)
 
 # 4. Test for differences in trajectory shape/progression
 # Kolmogorov-Smirnov test (tests if distributions are different)
 sex_levels <- levels(plot_df_clean$sex)
 ks_result <- ks.test(
   plot_df_clean$pseudotime[plot_df_clean$sex == sex_levels[1]],
   plot_df_clean$pseudotime[plot_df_clean$sex == sex_levels[2]]
 )
 print(paste("\nKS test p-value:", ks_result$p.value))
 
 # 5. Cell type composition along trajectory
 celltype_analysis <- plot_df_clean %>%
   mutate(pt_quartile = cut(pseudotime, 
                            breaks = quantile(pseudotime, probs = seq(0, 1, 0.25), na.rm = TRUE),
                            labels = c("Early", "Early-Mid", "Mid-Late", "Late"),
                            include.lowest = TRUE)) %>%
   filter(!is.na(pt_quartile)) %>%
   group_by(sex, pt_quartile, celltype) %>%
   summarise(n = n(), .groups = 'drop') %>%
   group_by(sex, pt_quartile) %>%
   mutate(prop = n / sum(n))
 
 p4 <- ggplot(celltype_analysis, aes(x = pt_quartile, y = prop, fill = celltype)) +
   geom_bar(stat = "identity", position = "stack") +
   facet_wrap(~sex) +
   labs(title = "Cell Type Composition Along Trajectory",
        x = "Trajectory Stage",
        y = "Proportion") +
   theme_classic() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
 print(p4)
 
 # ============================================
 # PART 2: PhenoPath Analysis (with enhanced debugging)
 # ============================================
 
 print("\n========================================")
 print("PART 2: PHENOPATH ANALYSIS")
 print("========================================\n")
 
 print("Running PhenoPath analysis...")
 
 # Get matching cells (handle potential mismatches)
 cells_to_keep <- intersect(rownames(plot_df_clean), colnames(so_subset))
 print(paste("Cells in both plot_df_clean and so_subset:", length(cells_to_keep)))
 
 # Subset Seurat object to matching cells
 so_subset_clean <- subset(so_subset, cells = cells_to_keep)
 
 # Prepare data for PhenoPath
 expr_matrix <- GetAssayData(so_subset_clean, slot = "data", assay = "RNA")
 
 # ENHANCED VALIDATION
 print("\n=== PhenoPath Input Validation ===")
 print(paste("Expression matrix dimensions:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "cells"))
 print(paste("Number of cells in so_subset_clean:", ncol(so_subset_clean)))
 
 # Check for NAs in expression matrix
 na_count <- sum(is.na(expr_matrix))
 print(paste("NAs in expression matrix:", na_count))
 
 # Convert sex to numeric (0/1)
 sex_numeric <- as.numeric(factor(so_subset_clean$sex)) - 1
 sex_levels_pheno <- levels(factor(so_subset_clean$sex))
 print(paste("\nSex encoding:", sex_levels_pheno[1], "= 0,", sex_levels_pheno[2], "= 1"))
 print(paste("Sex distribution:", sum(sex_numeric == 0), "vs", sum(sex_numeric == 1)))
 
 # Check for NAs in sex
 print(paste("NAs in sex_numeric:", sum(is.na(sex_numeric))))
 
 # Get pseudotime
 pseudotime_init <- so_subset_clean$pseudotime
 
 # CRITICAL: Check pseudotime
 print("\n=== Pseudotime Validation ===")
 print(paste("Length of pseudotime:", length(pseudotime_init)))
 print(paste("NAs in pseudotime:", sum(is.na(pseudotime_init))))
 print(paste("Range of pseudotime:", min(pseudotime_init, na.rm = TRUE), "to", max(pseudotime_init, na.rm = TRUE)))
 
 # If there are NAs in pseudotime, handle them
 if(any(is.na(pseudotime_init))) {
   print("WARNING: NAs found in pseudotime, filtering out these cells")
   
   valid_cells <- !is.na(pseudotime_init)
   
   so_subset_clean <- subset(so_subset_clean, cells = colnames(so_subset_clean)[valid_cells])
   expr_matrix <- GetAssayData(so_subset_clean, slot = "data", assay = "RNA")
   sex_numeric <- sex_numeric[valid_cells]
   pseudotime_init <- pseudotime_init[valid_cells]
   
   print(paste("After filtering: ", ncol(so_subset_clean), "cells remain"))
 }
 
 # Verify all vectors have same length
 print("\n=== Dimension Check ===")
 print(paste("Expression matrix cells:", ncol(expr_matrix)))
 print(paste("Sex vector length:", length(sex_numeric)))
 print(paste("Pseudotime vector length:", length(pseudotime_init)))
 
 if(!(ncol(expr_matrix) == length(sex_numeric) && length(sex_numeric) == length(pseudotime_init))) {
   stop("Dimension mismatch! All vectors must have the same length")
 }
 
 # Select genes for PhenoPath
 hvgs <- VariableFeatures(so_subset_clean)
 if(length(hvgs) > 2000) {
   genes_for_phenopath <- hvgs[1:2000]
 } else {
   genes_for_phenopath <- hvgs
 }
 
 print(paste("\nRunning PhenoPath on", length(genes_for_phenopath), "genes..."))
 
 # Filter genes with zero variance
 gene_vars <- apply(expr_matrix[genes_for_phenopath, ], 1, var)
 genes_with_var <- genes_for_phenopath[gene_vars > 0]
 print(paste("Genes with variance > 0:", length(genes_with_var)))
 
 if(length(genes_with_var) < 100) {
   stop("Too few genes with variance. Check your expression data.")
 }
 
 genes_for_phenopath <- genes_with_var
 
 # Prepare expression matrix for PhenoPath (cells x genes)
 expr_for_phenopath <- t(as.matrix(expr_matrix[genes_for_phenopath, ]))
 
 print("\n=== Final PhenoPath Input ===")
 print(paste("Expression matrix for PhenoPath:", nrow(expr_for_phenopath), "cells x", ncol(expr_for_phenopath), "genes"))
 print(paste("Sex covariate length:", length(sex_numeric)))
 print(paste("Pseudotime init length:", length(pseudotime_init)))
 
 # Check for any remaining NAs
 if(any(is.na(expr_for_phenopath))) {
   print("WARNING: NAs in expression matrix, replacing with 0")
   expr_for_phenopath[is.na(expr_for_phenopath)] <- 0
 }
 
 # Standardize pseudotime to [0, 1]
 pseudotime_scaled <- (pseudotime_init - min(pseudotime_init)) / (max(pseudotime_init) - min(pseudotime_init))
 print(paste("Scaled pseudotime range:", min(pseudotime_scaled), "to", max(pseudotime_scaled)))
 
 # Initialize flag
 phenopath_success <- FALSE
 
 # Try running PhenoPath with error handling
 tryCatch({
   print("\nAttempting to run PhenoPath...")
   
   phenopath_fit <- phenopath(
     exprs_obj = expr_for_phenopath,
     x = sex_numeric,
     elbo_tol = 1e-2,
     z_init = pseudotime_scaled,  # Use scaled pseudotime
     thin = 40,
     verbose = TRUE
   )
   
   print("PhenoPath completed successfully!")
   
   # Extract PhenoPath results
   phenopath_pseudotime <- phenopath_fit$z
   phenopath_beta <- phenopath_fit$beta
   phenopath_alpha <- phenopath_fit$alpha
   
   # Add to Seurat object
   so_subset_clean$phenopath_pseudotime <- phenopath_pseudotime
   
   # Create results data frame
   phenopath_results <- data.frame(
     gene = genes_for_phenopath,
     beta = phenopath_beta,
     alpha = phenopath_alpha,
     abs_beta = abs(phenopath_beta),
     abs_alpha = abs(phenopath_alpha)
   ) %>%
     arrange(desc(abs_beta))
   
   # Identify significant sex-interaction genes
   beta_threshold <- quantile(abs(phenopath_beta), 0.95)
   phenopath_results$sig_sex_effect <- abs(phenopath_results$beta) > beta_threshold
   
   print(paste("\nPhenoPath identified", sum(phenopath_results$sig_sex_effect), 
               "genes with strong sex-specific effects (top 5%)"))
   
   # Save PhenoPath results
   write.csv(phenopath_results, 
             "phenopath_sex_interaction_genes.csv", 
             row.names = FALSE)
   
   # Top sex-interaction genes from PhenoPath
   top_phenopath_genes <- head(phenopath_results, 20)
   print("\nTop 20 genes with sex-interaction effects (PhenoPath):")
   print(top_phenopath_genes[, c("gene", "beta", "alpha")])
   
   # Set flag that PhenoPath succeeded
   phenopath_success <- TRUE
   
 }, error = function(e) {
   print("ERROR in PhenoPath:")
   print(e)
   print("\nSkipping PhenoPath analysis and continuing with tradeSeq only...")
   
   # Create empty results so the rest of the script can continue
   phenopath_results <<- data.frame(
     gene = character(0),
     beta = numeric(0),
     alpha = numeric(0),
     abs_beta = numeric(0),
     abs_alpha = numeric(0),
     sig_sex_effect = logical(0)
   )
   
   phenopath_success <<- FALSE
 })
 
 # ============================================
 # PART 2B: Compare PhenoPath vs Slingshot (only if PhenoPath succeeded)
 # ============================================
 
 if(phenopath_success) {
   print("\n=== Comparing PhenoPath vs Slingshot ===")
   
   comparison_df <- data.frame(
     slingshot_pt = so_subset_clean$pseudotime,
     phenopath_pt = phenopath_pseudotime,
     sex = so_subset_clean$sex,
     celltype = so_subset_clean$KPMP_celltype,
     UMAP_1 = Embeddings(so_subset_clean, reduction = "umap")[, 1],
     UMAP_2 = Embeddings(so_subset_clean, reduction = "umap")[, 2]
   ) %>%
     filter(!is.na(slingshot_pt) & !is.na(phenopath_pt))
   
   # Correlation between methods
   cor_overall <- cor(comparison_df$slingshot_pt, comparison_df$phenopath_pt, 
                      use = "complete.obs")
   print(paste("Overall correlation between Slingshot and PhenoPath:", 
               round(cor_overall, 3)))
   
   # Correlation by sex
   cor_by_sex <- comparison_df %>%
     group_by(sex) %>%
     summarise(correlation = cor(slingshot_pt, phenopath_pt, use = "complete.obs"))
   print("Correlation by sex:")
   print(cor_by_sex)
   
   # Visualize comparison
   p5 <- ggplot(comparison_df, aes(x = slingshot_pt, y = phenopath_pt, color = sex)) +
     geom_point(alpha = 0.5, size = 0.8) +
     geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
     geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
     labs(title = "Slingshot vs PhenoPath Pseudotime",
          subtitle = paste("Overall r =", round(cor_overall, 3)),
          x = "Slingshot Pseudotime",
          y = "PhenoPath Pseudotime") +
     theme_classic() +
     facet_wrap(~sex)
   
   print(p5)
   
   # Visualize PhenoPath pseudotime on UMAP
   p6 <- ggplot(comparison_df, aes(x = UMAP_1, y = UMAP_2, color = phenopath_pt)) +
     geom_point(size = 0.5, alpha = 0.6) +
     facet_wrap(~sex) +
     scale_color_viridis(option = "magma") +
     labs(title = "PhenoPath Pseudotime by Sex",
          color = "PhenoPath\nPseudotime") +
     theme_classic()
   
   p7 <- ggplot(comparison_df, aes(x = UMAP_1, y = UMAP_2, color = slingshot_pt)) +
     geom_point(size = 0.5, alpha = 0.6) +
     facet_wrap(~sex) +
     scale_color_viridis(option = "magma") +
     labs(title = "Slingshot Pseudotime by Sex",
          color = "Slingshot\nPseudotime") +
     theme_classic()
   
   print(p6 / p7)
   
   # Test if PhenoPath pseudotime differs by sex more than Slingshot
   sling_sex_diff <- comparison_df %>%
     group_by(sex) %>%
     summarise(mean_pt = mean(slingshot_pt, na.rm = TRUE)) %>%
     pull(mean_pt) %>%
     diff() %>%
     abs()
   
   pheno_sex_diff <- comparison_df %>%
     group_by(sex) %>%
     summarise(mean_pt = mean(phenopath_pt, na.rm = TRUE)) %>%
     pull(mean_pt) %>%
     diff() %>%
     abs()
   
   print(paste("\nAbsolute difference in mean pseudotime between sexes:"))
   print(paste("  Slingshot:", round(sling_sex_diff, 4)))
   print(paste("  PhenoPath:", round(pheno_sex_diff, 4)))
   
   # Visualize top genes
   if(nrow(phenopath_results) > 0) {
     top_genes_phenopath <- head(phenopath_results$gene, 12)
     
     plots_phenopath <- lapply(top_genes_phenopath, function(gene) {
       gene_expr <- expr_matrix[gene, ]
       
       plot_data <- data.frame(
         expression = gene_expr,
         pseudotime = phenopath_pseudotime,
         sex = so_subset_clean$sex
       )
       
       ggplot(plot_data, aes(x = pseudotime, y = expression, color = sex)) +
         geom_point(alpha = 0.3, size = 0.5) +
         geom_smooth(method = "loess", se = TRUE) +
         labs(title = gene,
              subtitle = paste("β =", round(phenopath_results$beta[phenopath_results$gene == gene], 3)),
              x = "PhenoPath Pseudotime",
              y = "Expression") +
         theme_classic() +
         theme(legend.position = "bottom")
     })
     
     print(wrap_plots(plots_phenopath, ncol = 3))
   }
 } else {
   print("\nPhenoPath analysis was skipped due to errors. Continuing with tradeSeq only...")
 }
 
 # ============================================
 # PART 3: TradeSeq Analysis
 # ============================================
 
 print("\n========================================")
 print("PART 3: TRADESEQ ANALYSIS")
 print("========================================\n")
 
 # Prepare data for tradeSeq
 counts_matrix <- GetAssayData(so_subset_clean, slot = "counts", assay = "RNA")
 pseudotime_mat <- slingPseudotime(sling_res)
 cellweights_mat <- slingCurveWeights(sling_res)
 
 # Match cells
 cells_match <- intersect(colnames(counts_matrix), rownames(pseudotime_mat))
 counts_matrix <- counts_matrix[, cells_match]
 pseudotime_mat <- pseudotime_mat[cells_match, , drop = FALSE]
 cellweights_mat <- cellweights_mat[cells_match, , drop = FALSE]
 
 # Subset to genes tested
 genes_to_test <- genes_for_phenopath
 counts_matrix_subset <- counts_matrix[genes_to_test, ]
 
 print("Fitting GAM models with sex as covariate...")
 print("This may take several minutes...")
 
 # Fit GAM models with sex interaction
 set.seed(123)
 sce <- fitGAM(
   counts = counts_matrix_subset,
   pseudotime = pseudotime_mat,
   cellWeights = cellweights_mat,
   conditions = factor(so_subset_clean$sex[match(cells_match, colnames(so_subset_clean))]),
   nknots = 6,
   verbose = TRUE
 )
 
 print("Testing for sex-specific trajectory patterns...")
 
 # Test for condition (sex) effects
 sex_trajectory_test <- conditionTest(sce, global = TRUE, pairwise = TRUE)
 
 # Add gene names and sort by significance
 sex_trajectory_results <- sex_trajectory_test %>%
   as.data.frame() %>%
   mutate(gene = rownames(sex_trajectory_test),
          padj = p.adjust(pvalue, method = "BH"),
          significant = padj < 0.05) %>%
   arrange(pvalue)
 
 # Summary
 print(paste("\nTotal genes tested (tradeSeq):", nrow(sex_trajectory_results)))
 print(paste("Genes with sex-specific trajectory patterns (padj < 0.05):", 
             sum(sex_trajectory_results$significant)))
 print(paste("Genes with sex-specific trajectory patterns (padj < 0.01):", 
             sum(sex_trajectory_results$padj < 0.01)))
 
 # Top sex-differential genes
 top_sex_genes_tradeseq <- head(sex_trajectory_results, 20)
 print("\nTop 20 genes with sex-specific trajectory patterns (tradeSeq):")
 print(top_sex_genes_tradeseq[, c("gene", "pvalue", "padj")])
 
 # Save results
 write.csv(sex_trajectory_results, 
           "tradeseq_sex_specific_trajectory_genes.csv", 
           row.names = FALSE)
 
 # ============================================
 # PART 4: Compare PhenoPath and TradeSeq Results
 # ============================================
 
 if(phenopath_success && nrow(phenopath_results) > 0) {
   print("\n========================================")
   print("PART 4: METHOD COMPARISON")
   print("========================================\n")
   
   # Merge results from both methods
   combined_results <- phenopath_results %>%
     left_join(sex_trajectory_results, by = "gene", suffix = c("_phenopath", "_tradeseq"))
   
   # Compare rankings
   combined_results <- combined_results %>%
     mutate(
       phenopath_rank = rank(-abs_beta),
       tradeseq_rank = rank(pvalue)
     )
   
   # Correlation between methods
   rank_cor <- cor(combined_results$phenopath_rank, 
                   combined_results$tradeseq_rank, 
                   use = "complete.obs")
   print(paste("Rank correlation between PhenoPath and tradeSeq:", 
               round(rank_cor, 3)))
   
   # Scatter plot comparing methods
   p8 <- ggplot(combined_results, aes(x = abs_beta, y = -log10(pvalue))) +
     geom_point(alpha = 0.5, size = 1) +
     geom_vline(xintercept = beta_threshold, linetype = "dashed", color = "red") +
     geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
     labs(title = "Comparison: PhenoPath vs TradeSeq",
          subtitle = paste("Rank correlation =", round(rank_cor, 3)),
          x = "PhenoPath |β| (sex interaction)",
          y = "TradeSeq -log10(p-value)") +
     theme_classic()
   
   print(p8)
   
   # Identify genes significant in both methods
   sig_both <- combined_results %>%
     filter(sig_sex_effect == TRUE & significant == TRUE) %>%
     arrange(desc(abs_beta))
   
   print(paste("\nGenes identified by BOTH methods:", nrow(sig_both)))
   if(nrow(sig_both) > 0) {
     print("Top genes identified by both methods:")
     print(head(sig_both[, c("gene", "beta", "pvalue", "padj")], 10))
   }
   
   # Venn diagram of overlap
   venn.plot <- venn.diagram(
     x = list(
       PhenoPath = phenopath_results$gene[phenopath_results$sig_sex_effect],
       TradeSeq = sex_trajectory_results$gene[sex_trajectory_results$significant]
     ),
     category.names = c("PhenoPath", "TradeSeq"),
     filename = NULL,
     fill = c("#3498db", "#e74c3c"),
     alpha = 0.5
   )
   grid.draw(venn.plot)
   
   # Save combined results
   write.csv(combined_results, 
             "combined_phenopath_tradeseq_results.csv", 
             row.names = FALSE)
   
   # ============================================
   # PART 5: Visualize Top Genes from Both Methods
   # ============================================
   
   # Get top genes from each method
   top_phenopath_only <- head(phenopath_results$gene[phenopath_results$sig_sex_effect], 6)
   top_tradeseq_only <- head(sex_trajectory_results$gene[sex_trajectory_results$significant], 6)
   top_both_methods <- head(sig_both$gene, 6)
   
   # Combine for visualization
   genes_to_visualize <- unique(c(top_both_methods, top_phenopath_only[1:3], top_tradeseq_only[1:3]))
   
   # Plot with tradeSeq smoothers
   plots_comparison <- lapply(genes_to_visualize[1:min(9, length(genes_to_visualize))], function(gene) {
     plotSmoothers(sce, counts_matrix_subset, gene = gene, 
                   alpha = 1, border = TRUE) +
       ggtitle(paste(gene, 
                     "\nPhenoPath β:", round(combined_results$beta[combined_results$gene == gene], 3),
                     "| tradeSeq p:", signif(combined_results$pvalue[combined_results$gene == gene], 3)))
   })
   
   print(wrap_plots(plots_comparison, ncol = 3))
   
 } else {
   print("\n========================================")
   print("PART 4: Skipping method comparison (PhenoPath unavailable)")
   print("========================================\n")
 }
 # ============================================
 # PART 6: FOCUSED ANALYSIS - PT-S1/S2 SEX DIFFERENCES
 # ============================================
 
 print("\n========================================")
 print("PART 6: PT-S1/S2 SPECIFIC ANALYSIS")
 print("========================================\n")
 
 # Filter for PT-S1/S2 cells only
 pt_s1s2_df <- plot_df_clean %>%
   filter(celltype == "PT-S1/S2")
 
 print(paste("Total PT-S1/S2 cells:", nrow(pt_s1s2_df)))
 print("Distribution by sex:")
 print(table(pt_s1s2_df$sex))
 
 # 1. Statistical comparison of pseudotime in PT-S1/S2
 pt_s1s2_stats <- pt_s1s2_df %>%
   group_by(sex) %>%
   summarise(
     n_cells = n(),
     mean_pt = mean(pseudotime, na.rm = TRUE),
     median_pt = median(pseudotime, na.rm = TRUE),
     sd_pt = sd(pseudotime, na.rm = TRUE),
     se_pt = sd_pt / sqrt(n_cells),
     .groups = 'drop'
   )
 
 print("\nPT-S1/S2 Pseudotime Statistics by Sex:")
 print(pt_s1s2_stats)
 
 # Test pseudotime differences between sexes
 pt_s1s2_sex_counts <- table(pt_s1s2_df$sex)
 
 if(length(pt_s1s2_sex_counts) == 2 && all(pt_s1s2_sex_counts > 0)) {
   pt_s1s2_test <- wilcox.test(pseudotime ~ sex, data = pt_s1s2_df)
   print(paste("\nPT-S1/S2: Wilcoxon p-value =", signif(pt_s1s2_test$p.value, 3)))
   
   # T-test
   pt_s1s2_ttest <- t.test(pseudotime ~ sex, data = pt_s1s2_df)
   print(paste("PT-S1/S2: T-test p-value =", signif(pt_s1s2_ttest$p.value, 3)))
   
   # Effect size
   pt_s1s2_cohens <- cohen.d(pt_s1s2_df$pseudotime, pt_s1s2_df$sex)
   print(paste("PT-S1/S2: Cohen's d =", round(pt_s1s2_cohens$estimate, 3)))
 } else {
   print("\nPT-S1/S2: Cannot perform test - need both sexes represented")
   print(paste("Sex distribution:", paste(names(pt_s1s2_sex_counts), pt_s1s2_sex_counts, collapse = ", ")))
   pt_s1s2_test <- list(p.value = NA)
   pt_s1s2_ttest <- list(p.value = NA)
   pt_s1s2_cohens <- list(estimate = NA)
 }
 
 # 2. Visualizations for PT-S1/S2
 
 # Density plots by sex
 p_pt_density <- ggplot(pt_s1s2_df, aes(x = pseudotime, fill = sex)) +
   geom_density(alpha = 0.6) +
   labs(title = "Pseudotime Distribution in PT-S1/S2",
        subtitle = "Comparing Males vs Females",
        x = "Pseudotime",
        y = "Density") +
   theme_classic() +
   theme(legend.position = "bottom")
 
 print(p_pt_density)
 
 # Box plots with statistics
 p_pt_box <- ggplot(pt_s1s2_df, aes(x = sex, y = pseudotime, fill = sex)) +
   geom_violin(alpha = 0.5) +
   geom_boxplot(width = 0.3, alpha = 0.7, outlier.alpha = 0.3) +
   geom_jitter(width = 0.1, alpha = 0.1, size = 0.5) +
   stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
   labs(title = "Pseudotime in PT-S1/S2 by Sex",
        subtitle = "Diamond = mean, Box = median ± IQR",
        x = "Sex",
        y = "Pseudotime") +
   theme_classic() +
   theme(legend.position = "none")
 
 print(p_pt_box)
 
 # UMAP visualization for PT-S1/S2
 p_pt_umap <- ggplot(pt_s1s2_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
   geom_point(size = 1, alpha = 0.7) +
   facet_wrap(~sex) +
   scale_color_viridis(option = "magma") +
   labs(title = "PT-S1/S2 Pseudotime on UMAP",
        color = "Pseudotime") +
   theme_classic()
 
 print(p_pt_umap)
 
 # Ridge plot showing shift
 p_pt_ridge <- ggplot(pt_s1s2_df, aes(x = pseudotime, y = sex, fill = sex)) +
   geom_density_ridges(alpha = 0.7, scale = 1.5) +
   labs(title = "Pseudotime Distribution: PT-S1/S2 by Sex",
        x = "Pseudotime",
        y = "Sex") +
   theme_ridges() +
   theme(legend.position = "bottom")
 
 print(p_pt_ridge)
 
 # 3. Identify PT-S1/S2 specific sex-differential genes
 
 # Subset Seurat object to PT-S1/S2 cells
 pt_cells <- rownames(pt_s1s2_df)
 pt_cells_in_seurat <- intersect(pt_cells, colnames(so_subset_clean))
 so_pt <- subset(so_subset_clean, cells = pt_cells_in_seurat)
 
 print(paste("\nPT-S1/S2 subset contains", ncol(so_pt), "cells"))
 
 # Check if we have enough cells from both sexes for tradeSeq
 pt_sex_table <- table(so_pt$sex)
 print("Sex distribution in PT-S1/S2 subset:")
 print(pt_sex_table)
 
 if(length(pt_sex_table) == 2 && all(pt_sex_table >= 10)) {
   # Run tradeSeq specifically on PT-S1/S2 cells
   counts_pt <- GetAssayData(so_pt, slot = "counts", assay = "RNA")
   pseudotime_pt <- slingPseudotime(sling_res)[pt_cells_in_seurat, , drop = FALSE]
   cellweights_pt <- slingCurveWeights(sling_res)[pt_cells_in_seurat, , drop = FALSE]
   
   # Use HVGs or subset
   hvgs_pt <- VariableFeatures(so_pt)
   if(length(hvgs_pt) > 1500) {
     genes_pt <- hvgs_pt[1:1500]
   } else {
     genes_pt <- hvgs_pt
   }
   
   counts_pt_subset <- counts_pt[genes_pt, ]
   
   print("Fitting GAM models for PT-S1/S2 cells...")
   
   set.seed(456)
   sce_pt <- fitGAM(
     counts = counts_pt_subset,
     pseudotime = pseudotime_pt,
     cellWeights = cellweights_pt,
     conditions = factor(so_pt$sex),
     nknots = 6,
     verbose = TRUE
   )
   
   print("Testing for sex-specific patterns in PT-S1/S2...")
   sex_test_pt <- conditionTest(sce_pt, global = TRUE, pairwise = TRUE)
   
   sex_results_pt <- sex_test_pt %>%
     as.data.frame() %>%
     mutate(
       gene = rownames(sex_test_pt),
       padj = p.adjust(pvalue, method = "BH"),
       significant = padj < 0.05
     ) %>%
     arrange(pvalue)
   
   print(paste("\nPT-S1/S2: Genes with sex-specific patterns (padj < 0.05):", 
               sum(sex_results_pt$significant)))
   print(paste("PT-S1/S2: Genes with sex-specific patterns (padj < 0.01):", 
               sum(sex_results_pt$padj < 0.01)))
   
   # Save results
   write.csv(sex_results_pt, 
             "pt_s1s2_sex_specific_genes.csv", 
             row.names = FALSE)
   
   print("\nTop 20 PT-S1/S2 sex-differential genes:")
   print(head(sex_results_pt[, c("gene", "pvalue", "padj")], 20))
   
   # 4. Visualize top PT-S1/S2 sex-differential genes
   top_pt_genes <- head(sex_results_pt$gene, 12)
   
   plots_pt_genes <- lapply(top_pt_genes, function(gene) {
     plotSmoothers(sce_pt, counts_pt_subset, gene = gene, 
                   alpha = 1, border = TRUE) +
       ggtitle(paste(gene, "\nPT-S1/S2",
                     "(padj =", signif(sex_results_pt$padj[sex_results_pt$gene == gene], 3), ")"))
   })
   
   print(wrap_plots(plots_pt_genes, ncol = 3))
   
   # 5. Compare expression patterns between sexes
   expr_data <- GetAssayData(so_pt, slot = "data", assay = "RNA")
   
   pt_expr_summary <- expand.grid(
     gene = top_pt_genes[1:6],
     sex = levels(so_pt$sex),
     stringsAsFactors = FALSE
   )
   
   pt_expr_summary$mean_expr <- sapply(1:nrow(pt_expr_summary), function(i) {
     gene <- pt_expr_summary$gene[i]
     sx <- pt_expr_summary$sex[i]
     
     cells_subset <- colnames(so_pt)[so_pt$sex == sx]
     
     if(length(cells_subset) > 0) {
       mean(expr_data[gene, cells_subset])
     } else {
       NA
     }
   })
   
   # Plot expression patterns
   p_expr_pattern <- ggplot(pt_expr_summary, 
                            aes(x = sex, y = mean_expr, fill = sex)) +
     geom_bar(stat = "identity") +
     facet_wrap(~gene, scales = "free_y", ncol = 3) +
     labs(title = "Mean Expression in PT-S1/S2 by Sex",
          x = "Sex",
          y = "Mean Expression") +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1),
           legend.position = "none")
   
   print(p_expr_pattern)
   
   tradeseq_pt_success <- TRUE
   
 } else {
   print("\nSkipping tradeSeq for PT-S1/S2: insufficient cells in one or both sexes (need at least 10 per sex)")
   sex_results_pt <- data.frame(
     gene = character(0),
     pvalue = numeric(0),
     padj = numeric(0),
     significant = logical(0)
   )
   tradeseq_pt_success <- FALSE
 }
 
 # 6. Test if sex difference in pseudotime is specific to PT-S1/S2
 other_celltypes_df <- plot_df_clean %>%
   filter(celltype != "PT-S1/S2")
 
 if(nrow(other_celltypes_df) > 0 && nrow(pt_s1s2_df) > 0) {
   # Check if both have both sexes
   other_sex_counts <- table(other_celltypes_df$sex)
   
   if(length(other_sex_counts) == 2 && all(other_sex_counts > 0) && 
      length(pt_s1s2_sex_counts) == 2 && all(pt_s1s2_sex_counts > 0)) {
     
     other_test <- wilcox.test(pseudotime ~ sex, data = other_celltypes_df)
     
     print("\n\nComparison with other cell types:")
     print(paste("PT-S1/S2 p-value:", signif(pt_s1s2_test$p.value, 3)))
     print(paste("Other cell types p-value:", signif(other_test$p.value, 3)))
     
     # Visualization
     comparison_data <- bind_rows(
       pt_s1s2_df %>% mutate(group = "PT-S1/S2"),
       other_celltypes_df %>% mutate(group = "Other")
     )
     
     p_comparison <- ggplot(comparison_data, aes(x = sex, y = pseudotime, fill = sex)) +
       geom_violin(alpha = 0.5) +
       geom_boxplot(width = 0.3, alpha = 0.7) +
       facet_wrap(~group) +
       labs(title = "Sex Difference in Pseudotime: PT-S1/S2 vs Other Cells",
            x = "Sex",
            y = "Pseudotime") +
       theme_classic() +
       theme(legend.position = "none")
     
     print(p_comparison)
   }
 }
 
 # 7. Functional enrichment of PT-S1/S2 sex-differential genes
 if(tradeseq_pt_success && exists("sex_results_pt")) {
   sig_pt_genes <- sex_results_pt %>%
     filter(significant) %>%
     pull(gene)
   
   if(length(sig_pt_genes) > 10) {
     print("\nRunning GO enrichment for PT-S1/S2 genes...")
     
     tryCatch({
       gene_entrez_pt <- bitr(sig_pt_genes, 
                              fromType = "SYMBOL",
                              toType = "ENTREZID",
                              OrgDb = org.Hs.eg.db)
       
       go_results_pt <- enrichGO(
         gene = gene_entrez_pt$ENTREZID,
         OrgDb = org.Hs.eg.db,
         ont = "BP",
         pAdjustMethod = "BH",
         pvalueCutoff = 0.05,
         readable = TRUE
       )
       
       if(nrow(go_results_pt) > 0) {
         print("GO Enrichment for PT-S1/S2 sex-differential genes:")
         print(dotplot(go_results_pt, showCategory = 20) +
                 ggtitle("GO Enrichment: PT-S1/S2 Sex-Differential Genes"))
       }
       
       # KEGG pathway enrichment
       kegg_results_pt <- enrichKEGG(
         gene = gene_entrez_pt$ENTREZID,
         organism = "hsa",
         pAdjustMethod = "BH",
         pvalueCutoff = 0.05
       )
       
       if(nrow(kegg_results_pt) > 0) {
         print("KEGG Pathway Enrichment for PT-S1/S2 sex-differential genes:")
         print(dotplot(kegg_results_pt, showCategory = 15) +
                 ggtitle("KEGG Pathways: PT-S1/S2 Sex-Differential Genes"))
       }
     }, error = function(e) {
       print("Error in enrichment for PT-S1/S2 genes:")
       print(e)
     })
   }
 }
 
 # ============================================
 # PART 7: Functional Enrichment Analysis
 # ============================================
 
 print("\n========================================")
 print("PART 7: FUNCTIONAL ENRICHMENT")
 print("========================================\n")
 
 # Get significant genes from different analyses
 if(phenopath_success && nrow(phenopath_results) > 0) {
   sig_phenopath_genes <- phenopath_results %>%
     filter(sig_sex_effect) %>%
     pull(gene)
 } else {
   sig_phenopath_genes <- character(0)
 }
 
 sig_tradeseq_genes <- sex_trajectory_results %>%
   filter(significant) %>%
   pull(gene)
 
 sig_pt_genes <- sex_results_pt %>%
   filter(significant) %>%
   pull(gene)
 
 # Enrichment for PhenoPath genes
 if(length(sig_phenopath_genes) > 10) {
   print("Running GO enrichment for PhenoPath genes...")
   
   tryCatch({
     gene_entrez_pheno <- bitr(sig_phenopath_genes, 
                               fromType = "SYMBOL",
                               toType = "ENTREZID",
                               OrgDb = org.Hs.eg.db)
     
     go_results_pheno <- enrichGO(
       gene = gene_entrez_pheno$ENTREZID,
       OrgDb = org.Hs.eg.db,
       ont = "BP",
       pAdjustMethod = "BH",
       pvalueCutoff = 0.05,
       readable = TRUE
     )
     
     if(nrow(go_results_pheno) > 0) {
       print("GO Enrichment for PhenoPath genes:")
       print(dotplot(go_results_pheno, showCategory = 15) +
               ggtitle("GO Enrichment: PhenoPath Sex-Interaction Genes"))
     }
   }, error = function(e) {
     print("Error in GO enrichment for PhenoPath genes:")
     print(e)
   })
 }
 
 # Enrichment for tradeSeq genes
 if(length(sig_tradeseq_genes) > 10) {
   print("Running GO enrichment for tradeSeq genes...")
   
   tryCatch({
     gene_entrez_trade <- bitr(sig_tradeseq_genes, 
                               fromType = "SYMBOL",
                               toType = "ENTREZID",
                               OrgDb = org.Hs.eg.db)
     
     go_results_trade <- enrichGO(
       gene = gene_entrez_trade$ENTREZID,
       OrgDb = org.Hs.eg.db,
       ont = "BP",
       pAdjustMethod = "BH",
       pvalueCutoff = 0.05,
       readable = TRUE
     )
     
     if(nrow(go_results_trade) > 0) {
       print("GO Enrichment for tradeSeq genes:")
       print(dotplot(go_results_trade, showCategory = 15) +
               ggtitle("GO Enrichment: TradeSeq Sex-Differential Genes"))
     }
   }, error = function(e) {
     print("Error in GO enrichment for tradeSeq genes:")
     print(e)
   })
 }
 
 # Enrichment for PT-S1/S2 genes
 if(length(sig_pt_genes) > 10) {
   print("Running GO enrichment for PT-S1/S2 genes...")
   
   tryCatch({
     gene_entrez_pt <- bitr(sig_pt_genes, 
                            fromType = "SYMBOL",
                            toType = "ENTREZID",
                            OrgDb = org.Hs.eg.db)
     
     go_results_pt <- enrichGO(
       gene = gene_entrez_pt$ENTREZID,
       OrgDb = org.Hs.eg.db,
       ont = "BP",
       pAdjustMethod = "BH",
       pvalueCutoff = 0.05,
       readable = TRUE
     )
     
     if(nrow(go_results_pt) > 0) {
       print("GO Enrichment for PT-S1/S2 sex-differential genes:")
       print(dotplot(go_results_pt, showCategory = 20) +
               ggtitle("GO Enrichment: PT-S1/S2 Sex-Differential Genes"))
     }
     
     # KEGG pathway enrichment
     kegg_results_pt <- enrichKEGG(
       gene = gene_entrez_pt$ENTREZID,
       organism = "hsa",
       pAdjustMethod = "BH",
       pvalueCutoff = 0.05
     )
     
     if(nrow(kegg_results_pt) > 0) {
       print("KEGG Pathway Enrichment for PT-S1/S2 sex-differential genes:")
       print(dotplot(kegg_results_pt, showCategory = 15) +
               ggtitle("KEGG Pathways: PT-S1/S2 Sex-Differential Genes"))
     }
   }, error = function(e) {
     print("Error in enrichment for PT-S1/S2 genes:")
     print(e)
   })
 }
 
 # ============================================
 # FINAL SUMMARY REPORT
 # ============================================
 
 cat("\n========================================\n")
 cat("COMPREHENSIVE ANALYSIS SUMMARY\n")
 cat("========================================\n\n")
 
 cat("1. OVERALL TRAJECTORY COMPARISON:\n")
 cat(paste("   - Wilcoxon p-value:", signif(wilcox_result$p.value, 3), "\n"))
 cat(paste("   - Cohen's d:", round(cohens_d$estimate, 3), "\n"))
 cat(paste("   - KS test p-value:", signif(ks_result$p.value, 3), "\n\n"))
 
 if(phenopath_success) {
   cat("2. PHENOPATH RESULTS:\n")
   cat(paste("   - Genes with strong sex effects:", sum(phenopath_results$sig_sex_effect), "\n"))
   if(exists("cor_overall")) {
     cat(paste("   - Correlation with Slingshot:", round(cor_overall, 3), "\n\n"))
   }
 } else {
   cat("2. PHENOPATH RESULTS:\n")
   cat("   - PhenoPath analysis failed or was skipped\n\n")
 }
 
 cat("3. TRADESEQ RESULTS (ALL CELLS):\n")
 cat(paste("   - Significant genes (padj < 0.05):", sum(sex_trajectory_results$significant), "\n"))
 cat(paste("   - Significant genes (padj < 0.01):", sum(sex_trajectory_results$padj < 0.01), "\n\n"))
 
 if(phenopath_success && exists("combined_results")) {
   cat("4. METHOD COMPARISON:\n")
   cat(paste("   - Genes identified by PhenoPath:", sum(phenopath_results$sig_sex_effect), "\n"))
   cat(paste("   - Genes identified by tradeSeq:", sum(sex_trajectory_results$significant), "\n"))
   if(exists("sig_both")) {
     cat(paste("   - Genes identified by BOTH:", nrow(sig_both), "\n"))
   }
   if(exists("rank_cor")) {
     cat(paste("   - Rank correlation:", round(rank_cor, 3), "\n\n"))
   }
 } else {
   cat("4. METHOD COMPARISON:\n")
   cat("   - Skipped (PhenoPath unavailable)\n\n")
 }
 
 cat("5. PT-S1/S2 SPECIFIC ANALYSIS:\n")
 cat(paste("   - PT-S1 Wilcoxon p-value:", signif(pt_s1_test$p.value, 3), "\n"))
 cat(paste("   - PT-S1 Cohen's d:", round(pt_s1_cohens$estimate, 3), "\n"))
 cat(paste("   - PT-S2 Wilcoxon p-value:", signif(pt_s2_test$p.value, 3), "\n"))
 cat(paste("   - PT-S2 Cohen's d:", round(pt_s2_cohens$estimate, 3), "\n"))
 cat(paste("   - PT-S1/S2 sex-differential genes:", sum(sex_results_pt$significant), "\n\n"))
 
 cat("OUTPUT FILES CREATED:\n")
 if(phenopath_success) {
   cat("   - phenopath_sex_interaction_genes.csv\n")
 }
 cat("   - tradeseq_sex_specific_trajectory_genes.csv\n")
 if(phenopath_success && exists("combined_results")) {
   cat("   - combined_phenopath_tradeseq_results.csv\n")
 }
 cat("   - pt_s1s2_sex_specific_genes.csv\n")
 cat("   - Multiple visualization plots\n")
 
 cat("\n========================================\n")
 cat("Analysis complete!\n")
 cat("========================================\n\n")
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ##### Cell type proportions

 library(ggplot2)
 library(dplyr)
 library(tidyr)
 
 # Function to create stacked barplots by sex for each cell type
 create_celltype_barplots <- function(so_subset, dir.results, celltype, 
                                      subtype_column = "KPMP_celltype",
                                      test_method = "fisher",  # "fisher" or "proportion"
                                      alpha = 0.01) {  # More stringent alpha level
   
   # Subset Seurat object based on celltype
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
   
   # Clean celltype name for file saving
   celltype2 <- str_replace_all(celltype, "/", "_")
   celltype2 <- str_replace_all(celltype2, "-", "_")
   
   cat("\n=== Processing", celltype, "===\n")
   
   # Extract metadata
   metadata <- so_celltype@meta.data
   
   # Check if sex and subtype columns exist
   if(!"sex" %in% colnames(metadata)) {
     stop("'sex' column not found in metadata")
   }
   if(!subtype_column %in% colnames(metadata)) {
     stop(paste0("'", subtype_column, "' column not found in metadata"))
   }
   
   # Remove NAs
   metadata <- metadata %>%
     filter(!is.na(sex) & !is.na(!!sym(subtype_column)))
   
   # Calculate counts and proportions
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
   
   # Perform chi-square test for overall differences
   contingency_table <- metadata %>%
     count(sex, !!sym(subtype_column)) %>%
     pivot_wider(names_from = sex, values_from = n, values_fill = 0) %>%
     column_to_rownames(subtype_column) %>%
     as.matrix()
   
   chi_test <- chisq.test(contingency_table)
   
   # Only proceed with pairwise tests if overall chi-square is significant
   subtypes <- unique(metadata[[subtype_column]])
   prop_tests <- list()
   
   if(chi_test$p.value < 0.05) {  # Only do pairwise if overall is significant
     
     for(subtype in subtypes) {
       # Create 2x2 contingency table for this subtype
       subtype_data <- metadata %>%
         mutate(is_subtype = !!sym(subtype_column) == subtype)
       
       # Count for each sex
       test_table <- table(subtype_data$sex, subtype_data$is_subtype)
       
       # Only test if we have enough observations (at least 5 in each cell)
       if(all(test_table >= 5)) {
         
         if(test_method == "fisher") {
           # Fisher's exact test - more conservative, better for small samples
           stat_test <- fisher.test(test_table)
           p_val <- stat_test$p.value
         } else {
           # Proportion test
           stat_test <- prop.test(test_table)
           p_val <- stat_test$p.value
         }
         
         # Calculate the actual proportions for each sex
         male_prop <- test_table["Male", "TRUE"] / sum(test_table["Male", ])
         female_prop <- test_table["Female", "TRUE"] / sum(test_table["Female", ])
         
         # Calculate effect size (absolute difference)
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
   
   # Combine proportion test results
   if(length(prop_tests) > 0) {
     prop_test_results <- bind_rows(prop_tests)
     
     # Apply Bonferroni correction (more conservative than FDR)
     prop_test_results$p_bonferroni <- p.adjust(prop_test_results$p_value, method = "bonferroni")
     prop_test_results$p_fdr <- p.adjust(prop_test_results$p_value, method = "BH")
     
     # Require both: statistical significance AND meaningful effect size (>5% difference)
     prop_test_results$significant_bonferroni <- 
       (prop_test_results$p_bonferroni < alpha) & (prop_test_results$abs_difference > 5)
     
     prop_test_results$significant_fdr <- 
       (prop_test_results$p_fdr < alpha) & (prop_test_results$abs_difference > 5)
     
   } else {
     prop_test_results <- data.frame()
   }
   
   # Save statistical results
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
   
   # Add significance annotation to plot data (using Bonferroni)
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
   
   # Calculate midpoints for each subtype segment for asterisk placement
   cell_props <- cell_props %>%
     arrange(sex, desc(!!sym(subtype_column))) %>%
     group_by(sex) %>%
     mutate(
       cumsum_prop = cumsum(proportion),
       midpoint = cumsum_prop - proportion/2
     ) %>%
     ungroup()
   
   # Prepare asterisk data - only for Female bars (right bar)
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
   
   # Create color palette
   n_subtypes <- length(unique(cell_props[[subtype_column]]))
   colors <- scales::hue_pal()(n_subtypes)
   
   # Create the plot
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
   
   # Add asterisks for significant subtypes
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
   
   # Add caption with legend for asterisks
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
   
   # Save plots
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
 celltypes_vec <- c( 'PT', 'TAL', 'EC',  'DCTall', 'IC', 'intercalated')
 
 results_list <- list()
 for(celltype in celltypes_vec) {
   results_list[[celltype]] <- create_celltype_barplots(
     so_subset = so_subset, 
     dir.results = dir.results, 
     celltype = celltype,
     subtype_column = "KPMP_celltype",
     test_method = "fisher",  # Use "fisher" for more conservative, "proportion" for original
     alpha = 0.01  # More stringent threshold (0.01 instead of 0.05)
   )
   cat("\n")
 }
 
 
 
 
 