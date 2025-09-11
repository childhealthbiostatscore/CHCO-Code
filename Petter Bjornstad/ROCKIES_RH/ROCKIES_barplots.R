library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)







#barplots

lc_files <- list.files('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', pattern='csv')



t2d_files <- list.files('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2/', pattern = 'csv')


celltypes <- c('All', 'PT', 'PT-S1/S2', 'PT-S3', 'aPT', 'POD', 
               'TAL', 'C-TAL-1','C-TAL-2', 'dTAL', 'DCT', 'dDCT',
               'EC', 'EC-AEA', 'EC-AVR', 'EC-GC', 'EC-PTC', 
               "cDC",
               "cycT",
 #              "CD4+ T",
  #             "CD8+ T",
               "NK",
               "B",
               "MON",
               "MAC",
               "MC")


overall_summary <- data.frame(celltype = celltypes)
overall_summary$TCA_total <- 27
overall_summary$OxPhos_total <- 10
overall_summary$TCA_flipped <- NA
overall_summary$TCA_sig <- NA
overall_summary$OxPhos_flipped <- NA
overall_summary$OxPhos_sig <- NA


for(i in c(1:length(celltypes))){
  
  celltype <- celltypes[i]
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  tmp_lc <- lc_files[str_which(lc_files, pattern = paste0('cycle_', celltype2, '_cells'))]
  tmp_t2d <- t2d_files[str_which(t2d_files, pattern = paste0('cycle_', celltype2, '_cells'))]
  
  #TCA
  tmp_lc_tca <- tmp_lc[str_which(tmp_lc, pattern = 'TCA')]
  tmp_t2d_tca <- tmp_t2d[str_which(tmp_t2d, pattern = 'TCA')]
  
  tca_lc <- data.table::fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', tmp_lc_tca))
  tca_t2d <- data.table::fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2/', tmp_t2d_tca))
  
  tca_lc <- tca_lc %>% dplyr::select(gene, logFC_lc = logFC_groupType_2_Diabetes, 
                                     pvalue_lc =  any_of(c("p_groupType_2_Diabetes", "pvalue")))
  
  tca_t2d <- tca_t2d %>% dplyr::select(gene, logFC_t2d = logFC_epic_sglti2_1Yes, 
                                     pvalue_t2d =  any_of(c("p_epic_sglti2_1Yes", "pvalue")))
  
  tca_full <- tca_lc %>% left_join(tca_t2d)
  
  significant_traits <- tca_full %>% filter(pvalue_lc < 0.05 | pvalue_t2d < 0.05)
  opposite_1 <- tca_full %>% filter(logFC_lc > 0 & logFC_t2d < 0)
  opposite_2 <- tca_full %>% filter(logFC_lc < 0 & logFC_t2d > 0)
  
  tca_full_sig <- bind_rows(list(significant_traits, opposite_1, opposite_2)) %>% 
    filter(!duplicated(gene))
  
  overall_summary$TCA_flipped[i] <- tca_full_sig %>% filter((logFC_lc > 0 & logFC_t2d < 0) | (logFC_lc < 0 & logFC_t2d > 0)) %>% 
    nrow()
  overall_summary$TCA_sig[i] <- tca_full_sig %>% filter(pvalue_lc < 0.05 | pvalue_t2d < 0.05) %>% nrow()
  
  
  if(nrow(tca_full_sig) > 0){
    plot_df <- tca_full_sig %>% 
      gather(key = 'metric_condition', value = 'value', -gene) %>% 
      separate(metric_condition, into = c('metric', 'condition'), sep = '_') %>% 
      spread(key = metric, value = value) %>% mutate(Significance = ifelse(pvalue < 0.05, '*', ''))
    
    tmp_plot <- ggplot(plot_df, aes(x=gene, y=logFC, fill=condition))+
      geom_bar(stat='identity', position='dodge')+
      geom_text(aes(label = Significance), 
                position = position_dodge(width = 0.9),
                vjust = -0.5, 
                size = 4)+
      labs(x='Gene', y='LogFC', title = paste0('TCA Gene Comparison in ', celltype2, ' Cells'))+
      scale_fill_manual(name = 'Comparison',
                        labels = c('lc' = 'Lean Control vs. T2D (no SGLT2)',
                                   't2d' = 'T2D: no SGLT2 vs. SGLT2'),
                        values = c('lc' = 'darkturquoise', 
                                   't2d' = 'coral2'))+
      theme_classic()
    
    pdf(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/barplots/TCA_', celltype2, '_barplot.pdf'), 
        width = 12, height = 8)
    print(tmp_plot)
    dev.off()
    
    png(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/barplots/TCA_', celltype2, '_barplot.png'), 
        width = 1200, height = 800)
    print(tmp_plot)
    dev.off()
    
    
    
  }
  
  
  
  #OxPhos
  tmp_lc_oxphos <- tmp_lc[str_which(tmp_lc, pattern = 'PHOS_')]
  tmp_t2d_oxphos <- tmp_t2d[str_which(tmp_t2d, pattern = 'PHOS_')]
  
  oxphos_lc <- data.table::fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/', tmp_lc_oxphos))
  oxphos_t2d <- data.table::fread(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2/', tmp_t2d_oxphos))
  
  oxphos_lc <- oxphos_lc %>% dplyr::select(gene, logFC_lc = logFC_groupType_2_Diabetes, 
                                     pvalue_lc = any_of(c("p_groupType_2_Diabetes", "pvalue")))
  
  oxphos_t2d <- oxphos_t2d %>% dplyr::select(gene, logFC_t2d = logFC_epic_sglti2_1Yes, 
                                       pvalue_t2d = any_of(c("p_epic_sglti2_1Yes", "pvalue")))
  
  oxphos_full <- oxphos_lc %>% left_join(oxphos_t2d)
  
  significant_traits <- oxphos_full %>% filter(pvalue_lc < 0.05 | pvalue_t2d < 0.05)
  opposite_1 <- oxphos_full %>% filter(logFC_lc > 0 & logFC_t2d < 0)
  opposite_2 <- oxphos_full %>% filter(logFC_lc < 0 & logFC_t2d > 0)
  
  oxphos_full_sig <- bind_rows(list(significant_traits, opposite_1, opposite_2)) %>% 
    filter(!duplicated(gene))
  
  overall_summary$OxPhos_flipped[i] <- oxphos_full_sig %>% filter((logFC_lc > 0 & logFC_t2d < 0) | (logFC_lc < 0 & logFC_t2d > 0)) %>% 
    nrow()
  overall_summary$OxPhos_sig[i] <- oxphos_full_sig %>% filter(pvalue_lc < 0.05 | pvalue_t2d < 0.05) %>% nrow()
  
  if(nrow(oxphos_full_sig) > 0){
    plot_df <- oxphos_full_sig %>% 
      gather(key = 'metric_condition', value = 'value', -gene) %>% 
      separate(metric_condition, into = c('metric', 'condition'), sep = '_') %>% 
      spread(key = metric, value = value) %>% mutate(Significance = ifelse(pvalue < 0.05, '*', ''))
    
    tmp_plot <- ggplot(plot_df, aes(x=gene, y=logFC, fill=condition))+
      geom_bar(stat='identity', position='dodge')+
      geom_text(aes(label = Significance), 
                position = position_dodge(width = 0.9),
                vjust = -0.5, 
                size = 4)+
      labs(x='Gene', y='LogFC', title = paste0('OxPhos Gene Comparison in ', celltype2, ' Cells'))+
      scale_fill_manual(name = 'Comparison',
                        labels = c('lc' = 'Lean Control vs. T2D (no SGLT2)',
                                   't2d' = 'T2D: no SGLT2 vs. SGLT2'),
                        values = c('lc' = 'darkturquoise', 
                                   't2d' = 'coral2'))+
      theme_classic()
    
    pdf(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/barplots/oxphos_', celltype2, '_barplot.pdf'), 
        width = 12, height = 8)
    print(tmp_plot)
    dev.off()
    
    png(paste0('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/barplots/oxphos_', celltype2, '_barplot.png'), 
        width = 1200, height = 800)
    print(tmp_plot)
    dev.off()
    
    
    
    
  }
  
  
  
  
  print(paste0(celltype2, ' is done.'))
  
}



write.table(overall_summary, 'C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/barplots/Overall_SummaryofResults.txt', 
            row.names=F, quote=F, sep='\t')


overall_summary <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/barplots/Overall_SummaryofResults.txt')
library(gridExtra)


# Create summary data for the top panel (totals)
summary_data <- data.frame(
  Gene_Type = c("TCA", "OxPhos"),
  Total_Genes = c(27, 10)
)

# Prepare data for main plot (excluding "All" row for cell-type specific analysis)
plot_data <- overall_summary %>%
  pivot_longer(cols = c(TCA_flipped, TCA_sig, OxPhos_flipped, OxPhos_sig),
               names_to = "category", values_to = "count") %>%
  separate(category, into = c("gene_type", "status"), sep = "_") %>%
  mutate(
    gene_type = factor(gene_type, levels = c("TCA", "OxPhos")),
    status = factor(status, levels = c("flipped", "sig")),
    celltype = factor(celltype, levels = unique(celltype))
  )

main_plot <- ggplot(plot_data, aes(x = celltype, y = count, fill = interaction(gene_type, status))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  geom_hline(yintercept = 27, color = "#2E86C1", linetype = "solid", size = 1.2, alpha = 0.8) +
  geom_hline(yintercept = 10, color = "#E67E22", linetype = "solid", size = 1.2, alpha = 0.8) +
  annotate("text", x = 1, y = 28.2, label = "TCA Total (27)", color = "#2E86C1", 
           hjust = 0, fontface = "bold", size = 3.5) +
  annotate("text", x = 1, y = 11.2, label = "OxPhos Total (10)", color = "#E67E22", 
           hjust = 0, fontface = "bold", size = 3.5) +
  scale_fill_manual(
    name = "Category",
    values = c("TCA.flipped" = "#5DADE2", "TCA.sig" = "#F8C471",
               "OxPhos.flipped" = "#2E86C1", "OxPhos.sig" = "#E67E22"),
    labels = c("TCA Flipped", "TCA Significant", "OxPhos Flipped", "OxPhos Significant")
  ) +
  labs(
    title = "TCA and OxPhos Gene Analysis by Cell Type",
    x = "Cell Type",
    y = "Number of Genes"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))


# Alternative version: Save plots separately if needed
# ggsave("summary_plot.png", summary_plot, width = 8, height = 4, dpi = 300)
# ggsave("main_plot.png", main_plot, width = 12, height = 8, dpi = 300)

# Display the combined plot


pdf('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/barplots/OverallSummaryplot.pdf', width = 12, height=8)
print(main_plot)
dev.off()

png('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/barplots/OverallSummaryplot.png', width = 1200, height=800)
print(main_plot)
dev.off()



