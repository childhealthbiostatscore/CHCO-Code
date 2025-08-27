##ROCKIES Cell Distributions



combined <- data.table::fread('C:/Users/netio/Downloads/ROCKIES_scRNA_combined_LC_T2D.csv')

results.dir <- 'C:/Users/netio/Documents/UofW/Rockies/celltype_distributions/'



#Major cell types 

tmp <- combined %>% filter(celltype2 %in% c('PT', 'TAL', 'EC', 'POD'))
chisq.test(tmp$group2, tmp$celltype2)

png(paste0(results.dir, 'ROCKIES_majorcellypes_distributions.png'))
ggplot(tmp, aes(x=group2, fill=celltype2))+
  geom_bar(position='fill')+theme_classic()+labs(x='Condition Group', fill='Cell Type', title= 'Major Celltype Distributions')
dev.off()


#PT Cells
tmp <- combined %>% filter(KPMP_celltype %in% c('aPT', 'PT-S1/S2', 'PT-S3'))
chisq.test(tmp$group2, tmp$KPMP_celltype)

png(paste0(results.dir, 'ROCKIES_PTCells_distributions.png'))
ggplot(tmp, aes(x=group2, fill= KPMP_celltype))+
  geom_bar(position='fill')+theme_classic()+labs(x='Condition Group', fill='Cell Type', title= 'PT Cell Distributions')
dev.off()



#TAL Cells
tmp <- combined %>% filter(KPMP_celltype %in% c('C-TAL-1', 'C-TAL-2', 'dTAL'))
chisq.test(tmp$group2, tmp$KPMP_celltype)

png(paste0(results.dir, 'ROCKIES_TALCells_distributions.png'))
ggplot(tmp, aes(x=group2, fill= KPMP_celltype))+
  geom_bar(position='fill')+theme_classic()+labs(x='Condition Group', fill='Cell Type', title= 'TAL Cell Distributions')
dev.off()

#Endothelial Cells

tmp <- combined %>% filter(KPMP_celltype %in% c('EC-AEA', 'EC-AVR', 'EC-GC', 'EC-PTC'))
chisq.test(tmp$group2, tmp$KPMP_celltype)

png(paste0(results.dir, 'ROCKIES_EndothelialCells_distributions.png'))
ggplot(tmp, aes(x=group2, fill= KPMP_celltype))+
  geom_bar(position='fill')+theme_classic()+labs(x='Condition Group', fill='Cell Type', title= 'Endothelial Cell Distributions')
dev.off()


#Immune Cells

tmp <- combined %>% filter(KPMP_celltype %in% c("cDC",
                                                  "cycT",
                                                  "CD4+ T",
                                                  "CD8+ T",
                                                  "NK",
                                                  "B",
                                                  "MON",
                                                  "MAC",
                                                  "MC"))
chisq.test(tmp$group2, tmp$KPMP_celltype)

png(paste0(results.dir, 'ROCKIES_ImmuneCells_distributions.png'))
ggplot(tmp, aes(x=group2, fill= KPMP_celltype))+
  geom_bar(position='fill')+theme_classic()+labs(x='Condition Group', fill='Cell Type', title= 'Immune Cell Distributions')
dev.off()


























