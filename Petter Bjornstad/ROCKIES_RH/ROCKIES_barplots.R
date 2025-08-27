library(ggplot2)









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













