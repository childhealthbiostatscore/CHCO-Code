### PT Cell Subtype Analysis 



library(dplyr)
library(stringr)
library(ggplot2)
library(grid)
library(pheatmap)




#TCA Analysis 
lc_t2d_dir <- 'C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/'
t2d_sglt2_dir <- 'C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2'

lc_t2d_files <- list.files(lc_t2d_dir, pattern='TCA_cycle_a?PT')
t2d_sglt2_files <- list.files(t2d_sglt2_dir, pattern = 'TCA_cycle_a?PT')


for(){
  
}














#LC vs. T2D expression 



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


load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')

remove(so_kpmp_sc)

#dat_groups <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_GroupAssignments.txt')
#dat_groups <- dat_groups %>% filter(group2 %in% c('Lean Control', 'T2D-No SGLTi2'))

#so_subset <- subset(so_subset, record_id == dat_groups$record_id)
test <- so_subset@meta.data %>% dplyr::select(record_id, group) %>% filter(!duplicated(record_id))

#load('C:/Users/netio/Downloads/TCA_genes.txt')
#load('C:/Users/netio/Downloads/OxPhos_genes.txt')




genes_of_interest <- c('ACO1', 'ACO2', 'CS', 'DLAT', 'DLD', 'IDH1', 'IDH3A', 'MDH1', 'OGDH', 'PDHA1', 'PDHB', 'SDHB', 'SDHD', 'SUCLA2', 'SUCLG1', 'SUCLG2')

gene_data <- FetchData(so_subset, vars = genes_of_interest)
seurat_obj <- AddMetaData(so_subset, metadata = gene_data)


meta.data <- so_subset@meta.data

meta.data <- meta.data %>% dplyr::select(record_id, group, epic_sglti2_1, 
                                         celltype1, celltype2, KPMP_celltype, KPMP_celltype2) %>% 
  bind_cols(gene_data)


write.table(meta.data, 'C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/LC_T2D_metadata_withTCAgenes.txt', row.names=F, quote=F, 
            sep='\t')























#T2D SGLT2 expression 







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









load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')

dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2/'

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



counts_path <- round(GetAssayData(so_subset, layer = "counts")) # load counts and round


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



genes_of_interest <- c('ACO1', 'ACO2', 'CS', 'DLAT', 'DLD', 'IDH1', 'IDH3A', 'MDH1', 'OGDH', 'PDHA1', 'PDHB', 'SDHB', 'SDHD', 'SUCLA2', 'SUCLG1', 'SUCLG2')

gene_data <- FetchData(so_subset, vars = genes_of_interest)
seurat_obj <- AddMetaData(so_subset, metadata = gene_data)


meta.data <- so_subset@meta.data

meta.data <- meta.data %>% dplyr::select(record_id, group, epic_sglti2_1, 
                                         celltype1, celltype2, KPMP_celltype, KPMP_celltype2) %>% 
  bind_cols(gene_data)


write.table(meta.data, 'C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2/T2D_SGLT2_metadata_withTCAgenes.txt', row.names=F, quote=F, 
            sep='\t')







#Combined analysis 




t2d_meta <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/T2D_SGLT2/T2D_SGLT2_metadata_withTCAgenes.txt')
lc_t2d_meta <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/LC_T2D_metadata_withTCAgenes.txt') %>%
  anti_join(t2d_meta, by = 'record_id')

combined <- t2d_meta %>% bind_rows(lc_t2d_meta) %>% distinct() %>%
  mutate(group2 = ifelse(group == 'Lean_Control', 'LC', 
                          ifelse(epic_sglti2_1 == 'Yes', 'T2D-SGLT2', 'T2D-No SGLT2'))) %>%
  filter(KPMP_celltype %in% c('PT-S1/S2', 'PT-S3', 'aPT'))


combined_results <- combined %>% group_by(group2, KPMP_celltype) %>% 
  summarize(across(ACO1:SUCLG2, ~ mean(.x > 0) * 100, .names = 'pct_{.col}' ))

combined_results <- combined_results %>% mutate(axis_labels = paste0(KPMP_celltype, ': ', group2))

combined_results_long <- combined_results %>% 
  pivot_longer(
    cols = starts_with('pct'),
    names_to = 'gene', 
    values_to = 'percent_expressing') %>%
  mutate(gene = str_remove(gene, 'pct_'))














