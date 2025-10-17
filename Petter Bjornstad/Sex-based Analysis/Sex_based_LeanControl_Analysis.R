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

variable_names <- c('All', 'PT', 'TAL', 'EC', 'IC', 'DCTall')

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
celltypes_vec <- c('All', 'PT', 'TAL', 'EC', 'DCTall', 'IC')

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
   dplyr::select(record_id, group, sex, age, all_of(existing_cols))
 
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
 
 
 
 
 
 
 
 
 
 
 