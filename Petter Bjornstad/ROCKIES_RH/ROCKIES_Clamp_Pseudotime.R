########### Clamp and Pseudotime Analysis in ROCKIES 




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




### clamp

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

#date_of_screen
#screen_date

#dat <- harmonized_data %>% dplyr::select(-dob) %>% 
#  arrange(date_of_screen) %>% 
#  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
#                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
#                   .by = c(record_id, visit))



dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
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

dat2 <- dat %>% 
  dplyr::select(-avg_c_k2, -avg_c_f, -avg_m_k2, 
                -avg_m_f, -avg_c_k2)

dat2 <- dat2 %>% bind_cols(tmp_results)


dat2 <- dat2 %>% 
  dplyr::select(record_id, group, age, sex, bmi, m_i_p2_raw_lean, avg_c_k2, avg_c_f, avg_c_k2_f) %>% 
  filter(!is.na(m_i_p2_raw_lean)) %>% 
  filter(group != 'PKD') ##%>% 
 # filter(group == 'Lean Control')


cor(dat2 %>% dplyr::select(m_i_p2_raw_lean, avg_c_k2, avg_c_f, avg_c_k2_f), use = 'pairwise.complete.obs')


library(corrplot)
library(dplyr)

# Calculate correlation matrix
cor_matrix <- cor(dat2 %>% dplyr::select(m_i_p2_raw_lean, avg_c_k2, avg_c_f, avg_c_k2_f), 
                  use = 'pairwise.complete.obs')

# Rename rows and columns
rownames(cor_matrix) <- c("M/I", "Avg Cortical K2", "Avg Cortical F", "Avg Cortical K2/F")
colnames(cor_matrix) <- c("M/I", "Avg Cortical K2", "Avg Cortical F", "Avg Cortical K2/F")

# Extract just M/I row, excluding the M/I column (perfect correlation with itself)
cor_subset <- cor_matrix[1, 2:4, drop = FALSE]

# Save as PDF
pdf("C:/Users/netio/Documents/UofW/Rockies/correlation_plot_MI.pdf", width = 8, height = 3)
corrplot(cor_subset, 
         method = "color",           
         addCoef.col = "black",      
         number.cex = 1.2,           
         tl.col = "black",           
         tl.srt = 45,                
         col = colorRampPalette(c("blue", "white", "red"))(200),
         cl.pos = "b",               
         is.corr = TRUE)
dev.off()

# Save as PNG
png("C:/Users/netio/Documents/UofW/Rockies/correlation_plot_MI.png", width = 800, height = 300, res = 150)
corrplot(cor_subset, 
         method = "color",           
         addCoef.col = "black",      
         number.cex = 1.2,           
         tl.col = "black",           
         tl.srt = 45,                
         col = colorRampPalette(c("blue", "white", "red"))(200),
         cl.pos = "b",               
         is.corr = TRUE)
dev.off()

# Display in R
corrplot(cor_subset, 
         method = "color",           
         addCoef.col = "black",      
         number.cex = 1.2,           
         tl.col = "black",           
         tl.srt = 45,                
         col = colorRampPalette(c("blue", "white", "red"))(200),
         cl.pos = "b",               
         is.corr = TRUE)

















#########pseudotime 





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




#Get files 


load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')


so_subset <- so_kpmp_sc
remove(so_kpmp_sc)

test <- so_subset@meta.data

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


test$epic_sglti2_1 <- NULL

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

dat <- dat %>% 
  dplyr::select(-avg_c_k2, -avg_c_f)

dat <- dat %>% bind_cols(tmp_results)

dat <- dat %>% dplyr::select(record_id, epic_sglti2_1, avg_c_k2, avg_c_k2_f) %>% 
  filter(!is.na(epic_sglti2_1))



test <- test %>% left_join(dat, by='record_id')

so_subset@meta.data$epic_sglti2_1 <- test$epic_sglti2_1
so_subset@meta.data$avg_c_k2 <- test$avg_c_k2
so_subset@meta.data$avg_c_k2_f <- test$avg_c_k2_f
so_subset <- subset(so_subset, epic_sglti2_1 == 'No')


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


test <- so_subset@meta.data
test2 <- test %>% 
  dplyr::select(record_id, avg_c_k2, avg_c_k2_f, epic_sglti2_1) %>% 
  filter(!duplicated(record_id)) %>% 
mutate(
  avg_c_k2_f_binary = ifelse(avg_c_k2_f >= median(avg_c_k2_f, na.rm = TRUE), "high", "low"), 
  avg_c_k2_binary = ifelse(avg_c_k2 >= median(avg_c_k2, na.rm= TRUE), 'high', 'low')) %>% 
  dplyr::select(record_id, avg_c_k2_f_binary, avg_c_k2_binary)

test <- test %>% left_join(test2, by='record_id')

so_subset@meta.data$avg_c_k2_f_binary <- test$avg_c_k2_f_binary
so_subset@meta.data$avg_c_k2_binary <- test$avg_c_k2_binary

so_subset <- subset(so_subset, subset = avg_c_k2_f_binary %in% c('high', 'low'))


#Analysis
#PT Cells
so_subset <- subset(so_subset, subset = celltype2 == 'PT')
so_subset <- RunUMAP(so_subset, dims = 1:30)

sling_res <- slingshot(as.SingleCellExperiment(so_subset), clusterLabels = 'KPMP_celltype', 
                       start.clus = 'PT-S1', end.clus = 'aPT', reducedDim = 'UMAP')

so_subset$pseudotime <- slingPseudotime(sling_res)[,1]



dir.results <- 'C:/Users/netio/Documents/UofW/Rockies/pseudotime/'

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
  avg_c_k2_f_binary = so_subset$avg_c_k2_f_binary
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

# Plot B: Violin plot comparing pseudotime by Cortical K2/F level
p2 <- ggplot(plot_df_clean, aes(x = avg_c_k2_f_binary, y = pseudotime, fill = avg_c_k2_f_binary)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.1, size = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_classic() +
  labs(title = "Pseudotime Comparison by Cortical K2/F Level",
       x = "Cortical K2/F Level",
       y = "Pseudotime") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Plot C: Density plot by Cortical K2/F level
p3 <- ggplot(plot_df_clean, aes(x = pseudotime, fill = avg_c_k2_f_binary, color = avg_c_k2_f_binary)) +
  geom_density(alpha = 0.3, size = 1) +
  theme_classic() +
  labs(title = "Cell Density Along Pseudotime Trajectory by Cortical K2/F Level",
       x = "Pseudotime",
       y = "Density",
       fill = "Cortical K2/F",
       color = "Cortical K2/F") +
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

# Calculate total cell counts per cell type AND Cortical K2/F level
total_counts_per_celltype_k2f <- plot_df_clean %>%
  count(celltype, avg_c_k2_f_binary, name = "total_cells")

# Create data for all regions and calculate percentages BY Cortical K2/F level
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
  
  # Calculate counts in this region BY Cortical K2/F level
  region_counts <- region_cells %>%
    count(celltype, avg_c_k2_f_binary, name = "cells_in_region") %>%
    left_join(total_counts_per_celltype_k2f, by = c("celltype", "avg_c_k2_f_binary")) %>%
    mutate(
      percent_of_celltype = (cells_in_region / total_cells) * 100,
      region_number = region$number,
      region_name = region$name,
      region_start = region$start,
      region_end = region$end
    )
  
  # Add any missing cell type/Cortical K2/F combinations with 0%
  all_combinations <- expand.grid(
    celltype = unique(plot_df_clean$celltype),
    avg_c_k2_f_binary = unique(plot_df_clean$avg_c_k2_f_binary),
    stringsAsFactors = FALSE
  )
  
  existing_combinations <- region_counts %>% 
    select(celltype, avg_c_k2_f_binary) %>%
    distinct()
  
  missing_combinations <- all_combinations %>%
    anti_join(existing_combinations, by = c("celltype", "avg_c_k2_f_binary"))
  
  if(nrow(missing_combinations) > 0) {
    missing_data <- missing_combinations %>%
      left_join(total_counts_per_celltype_k2f, by = c("celltype", "avg_c_k2_f_binary")) %>%
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
cat("\n=== Percentage of Each Cell Type (by Cortical K2/F Level) in Each Region ===\n")
for(i in 1:3) {
  region_data <- all_region_data %>% filter(region_number == i)
  cat("\n", unique(region_data$region_name), "Region (Pseudotime", 
      round(unique(region_data$region_start), 2), "to", 
      round(unique(region_data$region_end), 2), "):\n")
  print(region_data %>% 
          select(celltype, avg_c_k2_f_binary, cells_in_region, total_cells, percent_of_celltype) %>%
          arrange(celltype, avg_c_k2_f_binary))
}

# Save the detailed composition data
write.csv(all_region_data, 
          paste0(dir.results, "Pseudotime_Regions_Celltype_Cortical_K2F_Percentage.csv"), 
          row.names = FALSE)

# Create bar plots showing percentage of each cell type in each region, SPLIT BY Cortical K2/F level
celltype_colors <- c("PT-S1/S2" = "#4DAF4A",  # Green
                     "PT-S3" = "#377EB8",      # Blue
                     "aPT" = "#E41A1C")        # Red

k2f_colors <- c("high" = "#E63946", "low" = "#457B9D")

bar_charts <- list()

for(i in 1:3) {
  region_data <- all_region_data %>% filter(region_number == i)
  region_name <- unique(region_data$region_name)
  region_start <- unique(region_data$region_start)
  region_end <- unique(region_data$region_end)
  
  bar_charts[[i]] <- ggplot(region_data, aes(x = celltype, y = percent_of_celltype, fill = avg_c_k2_f_binary)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.3) +
    scale_fill_manual(values = k2f_colors) +
    geom_text(aes(label = paste0(round(percent_of_celltype, 1), "%")),
              position = position_dodge(width = 0.9),
              vjust = -0.5, size = 3, fontface = "bold") +
    theme_classic() +
    labs(title = paste0(region_name, " Region\n(", round(region_start, 1), " - ", round(region_end, 1), ")"),
         x = "Cell Type",
         y = "% of Cell Type in Region",
         fill = "Cortical K2/F") +
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
ggsave(paste0(dir.results, "Complete_Pseudotime_Analysis_Cortical_K2F_with_Ridge.pdf"), 
       plot = combined_plot, width = 16, height = 16)
ggsave(paste0(dir.results, "Complete_Pseudotime_Analysis_Cortical_K2F_with_Ridge.png"), 
       plot = combined_plot, width = 16, height = 16, dpi = 300)


# Create a summary table showing key differences
summary_table <- all_region_data %>%
  select(region_name, celltype, avg_c_k2_f_binary, percent_of_celltype) %>%
  pivot_wider(names_from = avg_c_k2_f_binary, values_from = percent_of_celltype, values_fill = 0) %>%
  mutate(Difference = high - low) %>%
  arrange(region_name, celltype)

cat("\n=== Summary: High vs Low Cortical K2/F Differences by Region ===\n")
print(summary_table)

write.csv(summary_table, 
          paste0(dir.results, "Cortical_K2F_Differences_Summary.csv"), 
          row.names = FALSE)






#Pseudotime evaluation
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
  cortical_k2f = so_subset$avg_c_k2_f_binary  # Using avg_c_k2_f_binary
)



setwd('C:/Users/netio/Documents/UofW/Rockies/pseudotime/further_exploration/')

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
  cortical_k2f = so_subset$avg_c_k2_f_binary  # Using avg_c_k2_f_binary
)




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

print("Starting comprehensive Cortical K2/F-based trajectory analysis...")
print("========================================")

# Check for missing values in key variables
print("Checking for missing values:")
print(paste("Missing pseudotime:", sum(is.na(plot_df$pseudotime))))
print(paste("Missing Cortical K2/F:", sum(is.na(plot_df$cortical_k2f))))
print(paste("Missing UMAP_1:", sum(is.na(plot_df$UMAP_1))))
print(paste("Missing UMAP_2:", sum(is.na(plot_df$UMAP_2))))

# Remove rows with missing pseudotime or cortical_k2f
plot_df_clean <- plot_df %>%
  filter(!is.na(pseudotime) & !is.na(cortical_k2f) & !is.na(UMAP_1) & !is.na(UMAP_2))

print(paste("Original rows:", nrow(plot_df)))
print(paste("Clean rows:", nrow(plot_df_clean)))
print(paste("Removed rows:", nrow(plot_df) - nrow(plot_df_clean)))

# Check Cortical K2/F distribution
print("Cortical K2/F distribution:")
print(table(plot_df_clean$cortical_k2f))

# Make sure cortical_k2f is a factor with exactly 2 levels
plot_df_clean$cortical_k2f <- factor(plot_df_clean$cortical_k2f)
if(length(levels(plot_df_clean$cortical_k2f)) != 2) {
  stop("Cortical K2/F variable must have exactly 2 levels (high/low)")
}

# ============================================
# PART 1: Test for Overall Trajectory Differences
# ============================================

print("\n========================================")
print("PART 1: OVERALL TRAJECTORY ANALYSIS")
print("========================================\n")

# 1. Compare pseudotime distributions between Cortical K2/F levels
trajectory_stats <- plot_df_clean %>%
  group_by(cortical_k2f) %>%
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

print("Trajectory Statistics by Cortical K2/F:")
print(trajectory_stats)

# Save trajectory statistics
write.csv(trajectory_stats, "trajectory_stats_cortical_k2f.csv", row.names = FALSE)

# 2. Statistical tests for pseudotime differences
# Wilcoxon rank-sum test
wilcox_result <- wilcox.test(pseudotime ~ cortical_k2f, data = plot_df_clean)
print(paste("\nWilcoxon test p-value:", wilcox_result$p.value))

# T-test (if distributions are relatively normal)
t_result <- t.test(pseudotime ~ cortical_k2f, data = plot_df_clean)
print(paste("T-test p-value:", t_result$p.value))

# Effect size (Cohen's d)
cohens_d <- cohen.d(plot_df_clean$pseudotime, plot_df_clean$cortical_k2f)
print(paste("Cohen's d effect size:", cohens_d$estimate))

# 3. Visualize trajectory differences
p1 <- ggplot(plot_df_clean, aes(x = pseudotime, y = cortical_k2f, fill = cortical_k2f)) +
  geom_density_ridges(alpha = 0.7, scale = 0.9) +
  labs(title = "Pseudotime Distribution by Cortical K2/F Level",
       x = "Pseudotime",
       y = "Cortical K2/F") +
  theme_ridges() +
  theme(legend.position = "none")

p2 <- ggplot(plot_df_clean, aes(x = cortical_k2f, y = pseudotime, fill = cortical_k2f)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.5) +
  labs(title = "Pseudotime by Cortical K2/F Level",
       x = "Cortical K2/F",
       y = "Pseudotime") +
  theme_classic() +
  theme(legend.position = "none")

p3 <- ggplot(plot_df_clean, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.5, alpha = 0.6) +
  facet_wrap(~cortical_k2f) +
  scale_color_viridis(option = "magma") +
  labs(title = "Trajectory by Cortical K2/F Level on UMAP") +
  theme_classic()

# Combine and save plots
combined_p1_p3 <- (p1 / p2) | p3
print(combined_p1_p3)
ggsave("fig1_trajectory_comparison_cortical_k2f.pdf", combined_p1_p3, width = 14, height = 8)
ggsave("fig1_trajectory_comparison_cortical_k2f.png", combined_p1_p3, width = 14, height = 8, dpi = 300)

# 4. Test for differences in trajectory shape/progression
# Kolmogorov-Smirnov test (tests if distributions are different)
k2f_levels <- levels(plot_df_clean$cortical_k2f)
ks_result <- ks.test(
  plot_df_clean$pseudotime[plot_df_clean$cortical_k2f == k2f_levels[1]],
  plot_df_clean$pseudotime[plot_df_clean$cortical_k2f == k2f_levels[2]]
)
print(paste("\nKS test p-value:", ks_result$p.value))

# 5. Cell type composition along trajectory
celltype_analysis <- plot_df_clean %>%
  mutate(pt_quartile = cut(pseudotime, 
                           breaks = quantile(pseudotime, probs = seq(0, 1, 0.25), na.rm = TRUE),
                           labels = c("Early", "Early-Mid", "Mid-Late", "Late"),
                           include.lowest = TRUE)) %>%
  filter(!is.na(pt_quartile)) %>%
  group_by(cortical_k2f, pt_quartile, celltype) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(cortical_k2f, pt_quartile) %>%
  mutate(prop = n / sum(n))

# Save celltype analysis
write.csv(celltype_analysis, "celltype_composition_by_k2f.csv", row.names = FALSE)

p4 <- ggplot(celltype_analysis, aes(x = pt_quartile, y = prop, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~cortical_k2f) +
  labs(title = "Cell Type Composition Along Trajectory by Cortical K2/F Level",
       x = "Trajectory Stage",
       y = "Proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p4)
ggsave("fig2_celltype_composition_cortical_k2f.pdf", p4, width = 10, height = 6)
ggsave("fig2_celltype_composition_cortical_k2f.png", p4, width = 10, height = 6, dpi = 300)

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

# Convert Cortical K2/F to numeric (0/1)
k2f_numeric <- as.numeric(factor(so_subset_clean$avg_c_k2_f_binary)) - 1
k2f_levels_pheno <- levels(factor(so_subset_clean$avg_c_k2_f_binary))
print(paste("\nCortical K2/F encoding:", k2f_levels_pheno[1], "= 0,", k2f_levels_pheno[2], "= 1"))
print(paste("Cortical K2/F distribution:", sum(k2f_numeric == 0), "vs", sum(k2f_numeric == 1)))

# Check for NAs in Cortical K2/F
print(paste("NAs in k2f_numeric:", sum(is.na(k2f_numeric))))

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
  k2f_numeric <- k2f_numeric[valid_cells]
  pseudotime_init <- pseudotime_init[valid_cells]
  
  print(paste("After filtering: ", ncol(so_subset_clean), "cells remain"))
}

# Verify all vectors have same length
print("\n=== Dimension Check ===")
print(paste("Expression matrix cells:", ncol(expr_matrix)))
print(paste("Cortical K2/F vector length:", length(k2f_numeric)))
print(paste("Pseudotime vector length:", length(pseudotime_init)))

if(!(ncol(expr_matrix) == length(k2f_numeric) && length(k2f_numeric) == length(pseudotime_init))) {
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
print(paste("Cortical K2/F covariate length:", length(k2f_numeric)))
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
    x = k2f_numeric,
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
  
  # Identify significant Cortical K2/F-interaction genes
  beta_threshold <- quantile(abs(phenopath_beta), 0.95)
  phenopath_results$sig_k2f_effect <- abs(phenopath_results$beta) > beta_threshold
  
  print(paste("\nPhenoPath identified", sum(phenopath_results$sig_k2f_effect), 
              "genes with strong Cortical K2/F-specific effects (top 5%)"))
  
  # Save PhenoPath results
  write.csv(phenopath_results, 
            "phenopath_cortical_k2f_interaction_genes.csv", 
            row.names = FALSE)
  
  # Top Cortical K2/F-interaction genes from PhenoPath
  top_phenopath_genes <- head(phenopath_results, 20)
  print("\nTop 20 genes with Cortical K2/F-interaction effects (PhenoPath):")
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
    sig_k2f_effect = logical(0)
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
    cortical_k2f = so_subset_clean$avg_c_k2_f_binary,
    celltype = so_subset_clean$KPMP_celltype,
    UMAP_1 = Embeddings(so_subset_clean, reduction = "umap")[, 1],
    UMAP_2 = Embeddings(so_subset_clean, reduction = "umap")[, 2]
  ) %>%
    filter(!is.na(slingshot_pt) & !is.na(phenopath_pt))
  
  # Save comparison data
  write.csv(comparison_df, "phenopath_slingshot_comparison.csv", row.names = FALSE)
  
  # Correlation between methods
  cor_overall <- cor(comparison_df$slingshot_pt, comparison_df$phenopath_pt, 
                     use = "complete.obs")
  print(paste("Overall correlation between Slingshot and PhenoPath:", 
              round(cor_overall, 3)))
  
  # Correlation by Cortical K2/F level
  cor_by_k2f <- comparison_df %>%
    group_by(cortical_k2f) %>%
    summarise(correlation = cor(slingshot_pt, phenopath_pt, use = "complete.obs"))
  print("Correlation by Cortical K2/F level:")
  print(cor_by_k2f)
  write.csv(cor_by_k2f, "correlation_by_k2f_level.csv", row.names = FALSE)
  
  # Visualize comparison
  p5 <- ggplot(comparison_df, aes(x = slingshot_pt, y = phenopath_pt, color = cortical_k2f)) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    labs(title = "Slingshot vs PhenoPath Pseudotime",
         subtitle = paste("Overall r =", round(cor_overall, 3)),
         x = "Slingshot Pseudotime",
         y = "PhenoPath Pseudotime") +
    theme_classic() +
    facet_wrap(~cortical_k2f)
  
  print(p5)
  ggsave("fig3_slingshot_vs_phenopath_comparison.pdf", p5, width = 10, height = 5)
  ggsave("fig3_slingshot_vs_phenopath_comparison.png", p5, width = 10, height = 5, dpi = 300)
  
  # Visualize PhenoPath pseudotime on UMAP
  p6 <- ggplot(comparison_df, aes(x = UMAP_1, y = UMAP_2, color = phenopath_pt)) +
    geom_point(size = 0.5, alpha = 0.6) +
    facet_wrap(~cortical_k2f) +
    scale_color_viridis(option = "magma") +
    labs(title = "PhenoPath Pseudotime by Cortical K2/F Level",
         color = "PhenoPath\nPseudotime") +
    theme_classic()
  
  p7 <- ggplot(comparison_df, aes(x = UMAP_1, y = UMAP_2, color = slingshot_pt)) +
    geom_point(size = 0.5, alpha = 0.6) +
    facet_wrap(~cortical_k2f) +
    scale_color_viridis(option = "magma") +
    labs(title = "Slingshot Pseudotime by Cortical K2/F Level",
         color = "Slingshot\nPseudotime") +
    theme_classic()
  
  combined_p6_p7 <- p6 / p7
  print(combined_p6_p7)
  ggsave("fig4_pseudotime_umap_comparison.pdf", combined_p6_p7, width = 10, height = 10)
  ggsave("fig4_pseudotime_umap_comparison.png", combined_p6_p7, width = 10, height = 10, dpi = 300)
  
  # Test if PhenoPath pseudotime differs by Cortical K2/F more than Slingshot
  sling_k2f_diff <- comparison_df %>%
    group_by(cortical_k2f) %>%
    summarise(mean_pt = mean(slingshot_pt, na.rm = TRUE)) %>%
    pull(mean_pt) %>%
    diff() %>%
    abs()
  
  pheno_k2f_diff <- comparison_df %>%
    group_by(cortical_k2f) %>%
    summarise(mean_pt = mean(phenopath_pt, na.rm = TRUE)) %>%
    pull(mean_pt) %>%
    diff() %>%
    abs()
  
  print(paste("\nAbsolute difference in mean pseudotime between Cortical K2/F levels:"))
  print(paste("  Slingshot:", round(sling_k2f_diff, 4)))
  print(paste("  PhenoPath:", round(pheno_k2f_diff, 4)))
  
  # Visualize top genes
  if(nrow(phenopath_results) > 0) {
    top_genes_phenopath <- head(phenopath_results$gene, 12)
    
    plots_phenopath <- lapply(top_genes_phenopath, function(gene) {
      gene_expr <- expr_matrix[gene, ]
      
      plot_data <- data.frame(
        expression = gene_expr,
        pseudotime = phenopath_pseudotime,
        cortical_k2f = so_subset_clean$avg_c_k2_f_binary
      )
      
      ggplot(plot_data, aes(x = pseudotime, y = expression, color = cortical_k2f)) +
        geom_point(alpha = 0.3, size = 0.5) +
        geom_smooth(method = "loess", se = TRUE) +
        labs(title = gene,
             subtitle = paste("β =", round(phenopath_results$beta[phenopath_results$gene == gene], 3)),
             x = "PhenoPath Pseudotime",
             y = "Expression") +
        theme_classic() +
        theme(legend.position = "bottom")
    })
    
    combined_phenopath_genes <- wrap_plots(plots_phenopath, ncol = 3)
    print(combined_phenopath_genes)
    ggsave("fig5_top_phenopath_genes.pdf", combined_phenopath_genes, width = 15, height = 12)
    ggsave("fig5_top_phenopath_genes.png", combined_phenopath_genes, width = 15, height = 12, dpi = 300)
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

print("Fitting GAM models with Cortical K2/F as covariate...")
print("This may take several minutes...")

# Fit GAM models with Cortical K2/F interaction
set.seed(123)
sce <- fitGAM(
  counts = counts_matrix_subset,
  pseudotime = pseudotime_mat,
  cellWeights = cellweights_mat,
  conditions = factor(so_subset_clean$avg_c_k2_f_binary[match(cells_match, colnames(so_subset_clean))]),
  nknots = 6,
  verbose = TRUE
)

print("Testing for Cortical K2/F-specific trajectory patterns...")

# Test for condition (Cortical K2/F) effects
k2f_trajectory_test <- conditionTest(sce, global = TRUE, pairwise = TRUE)

# Add gene names and sort by significance
k2f_trajectory_results <- k2f_trajectory_test %>%
  as.data.frame() %>%
  mutate(gene = rownames(k2f_trajectory_test),
         padj = p.adjust(pvalue, method = "BH"),
         significant = padj < 0.05) %>%
  arrange(pvalue)

# Summary
print(paste("\nTotal genes tested (tradeSeq):", nrow(k2f_trajectory_results)))
print(paste("Genes with Cortical K2/F-specific trajectory patterns (padj < 0.05):", 
            sum(k2f_trajectory_results$significant)))
print(paste("Genes with Cortical K2/F-specific trajectory patterns (padj < 0.01):", 
            sum(k2f_trajectory_results$padj < 0.01)))

# Top Cortical K2/F-differential genes
top_k2f_genes_tradeseq <- head(k2f_trajectory_results, 20)
print("\nTop 20 genes with Cortical K2/F-specific trajectory patterns (tradeSeq):")
print(top_k2f_genes_tradeseq[, c("gene", "pvalue", "padj")])

# Save results
write.csv(k2f_trajectory_results, 
          "tradeseq_cortical_k2f_specific_trajectory_genes.csv", 
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
    left_join(k2f_trajectory_results, by = "gene", suffix = c("_phenopath", "_tradeseq"))
  
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
         x = "PhenoPath |β| (Cortical K2/F interaction)",
         y = "TradeSeq -log10(p-value)") +
    theme_classic()
  
  print(p8)
  ggsave("fig6_phenopath_vs_tradeseq_comparison.pdf", p8, width = 8, height = 6)
  ggsave("fig6_phenopath_vs_tradeseq_comparison.png", p8, width = 8, height = 6, dpi = 300)
  
  # Identify genes significant in both methods
  sig_both <- combined_results %>%
    filter(sig_k2f_effect == TRUE & significant == TRUE) %>%
    arrange(desc(abs_beta))
  
  print(paste("\nGenes identified by BOTH methods:", nrow(sig_both)))
  if(nrow(sig_both) > 0) {
    print("Top genes identified by both methods:")
    print(head(sig_both[, c("gene", "beta", "pvalue", "padj")], 10))
  }
  
  # Venn diagram of overlap
  venn.plot <- venn.diagram(
    x = list(
      PhenoPath = phenopath_results$gene[phenopath_results$sig_k2f_effect],
      TradeSeq = k2f_trajectory_results$gene[k2f_trajectory_results$significant]
    ),
    category.names = c("PhenoPath", "TradeSeq"),
    filename = NULL,
    fill = c("#3498db", "#e74c3c"),
    alpha = 0.5
  )
  
  pdf("fig7_venn_diagram_methods_overlap.pdf", width = 8, height = 8)
  grid.draw(venn.plot)
  dev.off()
  
  png("fig7_venn_diagram_methods_overlap.png", width = 8, height = 8, units = "in", res = 300)
  grid.draw(venn.plot)
  dev.off()
  
  grid.draw(venn.plot)
  
  # Save combined results
  write.csv(combined_results, 
            "combined_phenopath_tradeseq_cortical_k2f_results.csv", 
            row.names = FALSE)
  
  # ============================================
  # PART 5: Visualize Top Genes from Both Methods
  # ============================================
  
  # Get top genes from each method
  top_phenopath_only <- head(phenopath_results$gene[phenopath_results$sig_k2f_effect], 6)
  top_tradeseq_only <- head(k2f_trajectory_results$gene[k2f_trajectory_results$significant], 6)
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
  
  combined_comparison_plots <- wrap_plots(plots_comparison, ncol = 3)
  print(combined_comparison_plots)
  ggsave("fig8_top_genes_both_methods.pdf", combined_comparison_plots, width = 15, height = 12)
  ggsave("fig8_top_genes_both_methods.png", combined_comparison_plots, width = 15, height = 12, dpi = 300)
  
} else {
  print("\n========================================")
  print("PART 4: Skipping method comparison (PhenoPath unavailable)")
  print("========================================\n")
}

# ============================================
# PART 6: FOCUSED ANALYSIS - ALL PT CELLS
# ============================================

print("\n========================================")
print("PART 6: ALL PT CELLS SPECIFIC ANALYSIS")
print("========================================\n")

# Filter for all PT cells (PT-S1/S2, PT-S3, aPT)
pt_all_df <- plot_df_clean %>%
  filter(celltype %in% c("PT-S1/S2", "PT-S3", "aPT"))

print(paste("Total PT cells:", nrow(pt_all_df)))
print("Distribution by Cortical K2/F level:")
print(table(pt_all_df$cortical_k2f))
print("\nDistribution by cell type:")
print(table(pt_all_df$celltype))

# 1. Statistical comparison of pseudotime in all PT cells
pt_all_stats <- pt_all_df %>%
  group_by(cortical_k2f) %>%
  summarise(
    n_cells = n(),
    mean_pt = mean(pseudotime, na.rm = TRUE),
    median_pt = median(pseudotime, na.rm = TRUE),
    sd_pt = sd(pseudotime, na.rm = TRUE),
    se_pt = sd_pt / sqrt(n_cells),
    .groups = 'drop'
  )

print("\nAll PT Cells Pseudotime Statistics by Cortical K2/F:")
print(pt_all_stats)
write.csv(pt_all_stats, "pt_all_cells_stats_k2f.csv", row.names = FALSE)

# Test pseudotime differences between Cortical K2/F levels
pt_all_k2f_counts <- table(pt_all_df$cortical_k2f)

if(length(pt_all_k2f_counts) == 2 && all(pt_all_k2f_counts > 0)) {
  pt_all_test <- wilcox.test(pseudotime ~ cortical_k2f, data = pt_all_df)
  print(paste("\nAll PT Cells: Wilcoxon p-value =", signif(pt_all_test$p.value, 3)))
  
  # T-test
  pt_all_ttest <- t.test(pseudotime ~ cortical_k2f, data = pt_all_df)
  print(paste("All PT Cells: T-test p-value =", signif(pt_all_ttest$p.value, 3)))
  
  # Effect size
  pt_all_cohens <- cohen.d(pt_all_df$pseudotime, pt_all_df$cortical_k2f)
  print(paste("All PT Cells: Cohen's d =", round(pt_all_cohens$estimate, 3)))
} else {
  print("\nAll PT Cells: Cannot perform test - need both Cortical K2/F levels represented")
  print(paste("Cortical K2/F distribution:", paste(names(pt_all_k2f_counts), pt_all_k2f_counts, collapse = ", ")))
  pt_all_test <- list(p.value = NA)
  pt_all_ttest <- list(p.value = NA)
  pt_all_cohens <- list(estimate = NA)
}

# 2. Visualizations for all PT cells

# Density plots by Cortical K2/F level
p_pt_density <- ggplot(pt_all_df, aes(x = pseudotime, fill = cortical_k2f)) +
  geom_density(alpha = 0.6) +
  labs(title = "Pseudotime Distribution in All PT Cells",
       subtitle = "Comparing High vs Low Cortical K2/F",
       x = "Pseudotime",
       y = "Density") +
  theme_classic() +
  theme(legend.position = "bottom")

print(p_pt_density)
ggsave("fig9_pt_all_density_k2f.pdf", p_pt_density, width = 8, height = 6)
ggsave("fig9_pt_all_density_k2f.png", p_pt_density, width = 8, height = 6, dpi = 300)

# Box plots with statistics
p_pt_box <- ggplot(pt_all_df, aes(x = cortical_k2f, y = pseudotime, fill = cortical_k2f)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.3, alpha = 0.7, outlier.alpha = 0.3) +
  geom_jitter(width = 0.1, alpha = 0.1, size = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  labs(title = "Pseudotime in All PT Cells by Cortical K2/F",
       subtitle = "Diamond = mean, Box = median ± IQR",
       x = "Cortical K2/F",
       y = "Pseudotime") +
  theme_classic() +
  theme(legend.position = "none")

print(p_pt_box)
ggsave("fig10_pt_all_boxplot_k2f.pdf", p_pt_box, width = 6, height = 6)
ggsave("fig10_pt_all_boxplot_k2f.png", p_pt_box, width = 6, height = 6, dpi = 300)

# UMAP visualization for all PT cells
p_pt_umap <- ggplot(pt_all_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 1, alpha = 0.7) +
  facet_wrap(~cortical_k2f) +
  scale_color_viridis(option = "magma") +
  labs(title = "All PT Cells Pseudotime on UMAP",
       color = "Pseudotime") +
  theme_classic()

print(p_pt_umap)
ggsave("fig11_pt_all_umap_k2f.pdf", p_pt_umap, width = 10, height = 5)
ggsave("fig11_pt_all_umap_k2f.png", p_pt_umap, width = 10, height = 5, dpi = 300)

# Ridge plot showing shift
p_pt_ridge <- ggplot(pt_all_df, aes(x = pseudotime, y = cortical_k2f, fill = cortical_k2f)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5) +
  labs(title = "Pseudotime Distribution: All PT Cells by Cortical K2/F",
       x = "Pseudotime",
       y = "Cortical K2/F") +
  theme_ridges() +
  theme(legend.position = "bottom")

print(p_pt_ridge)
ggsave("fig12_pt_all_ridge_k2f.pdf", p_pt_ridge, width = 8, height = 6)
ggsave("fig12_pt_all_ridge_k2f.png", p_pt_ridge, width = 8, height = 6, dpi = 300)

# Faceted by cell type
p_pt_by_celltype <- ggplot(pt_all_df, aes(x = pseudotime, fill = cortical_k2f)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~celltype, ncol = 1) +
  labs(title = "Pseudotime Distribution by PT Cell Type and Cortical K2/F",
       x = "Pseudotime",
       y = "Density",
       fill = "Cortical K2/F") +
  theme_classic() +
  theme(legend.position = "bottom")

print(p_pt_by_celltype)
ggsave("fig13_pt_by_celltype_k2f.pdf", p_pt_by_celltype, width = 8, height = 10)
ggsave("fig13_pt_by_celltype_k2f.png", p_pt_by_celltype, width = 8, height = 10, dpi = 300)

# 3. Identify PT-specific Cortical K2/F-differential genes

# Subset Seurat object to all PT cells
pt_cells <- rownames(pt_all_df)
pt_cells_in_seurat <- intersect(pt_cells, colnames(so_subset_clean))
so_pt <- subset(so_subset_clean, cells = pt_cells_in_seurat)

print(paste("\nAll PT subset contains", ncol(so_pt), "cells"))

# Check if we have enough cells from both Cortical K2/F levels for tradeSeq
pt_k2f_table <- table(so_pt$avg_c_k2_f_binary)
print("Cortical K2/F distribution in all PT subset:")
print(pt_k2f_table)

if(length(pt_k2f_table) == 2 && all(pt_k2f_table >= 10)) {
  # Run tradeSeq specifically on all PT cells
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
  
  print("Fitting GAM models for all PT cells...")
  
  set.seed(456)
  sce_pt <- fitGAM(
    counts = counts_pt_subset,
    pseudotime = pseudotime_pt,
    cellWeights = cellweights_pt,
    conditions = factor(so_pt$avg_c_k2_f_binary),
    nknots = 6,
    verbose = TRUE
  )
  
  print("Testing for Cortical K2/F-specific patterns in all PT cells...")
  k2f_test_pt <- conditionTest(sce_pt, global = TRUE, pairwise = TRUE)
  
  k2f_results_pt <- k2f_test_pt %>%
    as.data.frame() %>%
    mutate(
      gene = rownames(k2f_test_pt),
      padj = p.adjust(pvalue, method = "BH"),
      significant = padj < 0.05
    ) %>%
    arrange(pvalue)
  
  print(paste("\nAll PT Cells: Genes with Cortical K2/F-specific patterns (padj < 0.05):", 
              sum(k2f_results_pt$significant)))
  print(paste("All PT Cells: Genes with Cortical K2/F-specific patterns (padj < 0.01):", 
              sum(k2f_results_pt$padj < 0.01)))
  
  # Save results
  write.csv(k2f_results_pt, 
            "pt_all_cells_cortical_k2f_specific_genes.csv", 
            row.names = FALSE)
  
  print("\nTop 20 All PT Cortical K2/F-differential genes:")
  print(head(k2f_results_pt[, c("gene", "pvalue", "padj")], 20))
  
  # 4. Visualize top all PT Cortical K2/F-differential genes
  top_pt_genes <- head(k2f_results_pt$gene, 12)
  
  plots_pt_genes <- lapply(top_pt_genes, function(gene) {
    plotSmoothers(sce_pt, counts_pt_subset, gene = gene, 
                  alpha = 1, border = TRUE) +
      ggtitle(paste(gene, "\nAll PT Cells",
                    "(padj =", signif(k2f_results_pt$padj[k2f_results_pt$gene == gene], 3), ")"))
  })
  
  combined_pt_genes <- wrap_plots(plots_pt_genes, ncol = 3)
  print(combined_pt_genes)
  ggsave("fig14_top_pt_all_k2f_genes.pdf", combined_pt_genes, width = 15, height = 12)
  ggsave("fig14_top_pt_all_k2f_genes.png", combined_pt_genes, width = 15, height = 12, dpi = 300)
  
  # 5. Compare expression patterns between Cortical K2/F levels
  expr_data <- GetAssayData(so_pt, slot = "data", assay = "RNA")
  
  pt_expr_summary <- expand.grid(
    gene = top_pt_genes[1:6],
    cortical_k2f = levels(factor(so_pt$avg_c_k2_f_binary)),
    stringsAsFactors = FALSE
  )
  
  pt_expr_summary$mean_expr <- sapply(1:nrow(pt_expr_summary), function(i) {
    gene <- pt_expr_summary$gene[i]
    k2f <- pt_expr_summary$cortical_k2f[i]
    
    cells_subset <- colnames(so_pt)[so_pt$avg_c_k2_f_binary == k2f]
    
    if(length(cells_subset) > 0) {
      mean(expr_data[gene, cells_subset])
    } else {
      NA
    }
  })
  
  # Save expression summary
  write.csv(pt_expr_summary, "pt_all_expr_summary_k2f.csv", row.names = FALSE)
  
  # Plot expression patterns
  p_expr_pattern <- ggplot(pt_expr_summary, 
                           aes(x = cortical_k2f, y = mean_expr, fill = cortical_k2f)) +
    geom_bar(stat = "identity") +
    facet_wrap(~gene, scales = "free_y", ncol = 3) +
    labs(title = "Mean Expression in All PT Cells by Cortical K2/F",
         x = "Cortical K2/F",
         y = "Mean Expression") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  print(p_expr_pattern)
  ggsave("fig15_pt_all_expr_pattern_k2f.pdf", p_expr_pattern, width = 10, height = 8)
  ggsave("fig15_pt_all_expr_pattern_k2f.png", p_expr_pattern, width = 10, height = 8, dpi = 300)
  
  tradeseq_pt_success <- TRUE
  
} else {
  print("\nSkipping tradeSeq for all PT cells: insufficient cells in one or both Cortical K2/F levels (need at least 10 per level)")
  k2f_results_pt <- data.frame(
    gene = character(0),
    pvalue = numeric(0),
    padj = numeric(0),
    significant = logical(0)
  )
  tradeseq_pt_success <- FALSE
}

# 6. Test if Cortical K2/F difference in pseudotime is specific to PT cells
other_celltypes_df <- plot_df_clean %>%
  filter(!celltype %in% c("PT-S1/S2", "PT-S3", "aPT"))

if(nrow(other_celltypes_df) > 0 && nrow(pt_all_df) > 0) {
  # Check if both have both Cortical K2/F levels
  other_k2f_counts <- table(other_celltypes_df$cortical_k2f)
  
  if(length(other_k2f_counts) == 2 && all(other_k2f_counts > 0) && 
     length(pt_all_k2f_counts) == 2 && all(pt_all_k2f_counts > 0)) {
    
    other_test <- wilcox.test(pseudotime ~ cortical_k2f, data = other_celltypes_df)
    
    print("\n\nComparison with non-PT cell types:")
    print(paste("All PT cells p-value:", signif(pt_all_test$p.value, 3)))
    print(paste("Non-PT cell types p-value:", signif(other_test$p.value, 3)))
    
    # Visualization
    comparison_data <- bind_rows(
      pt_all_df %>% mutate(group = "All PT Cells"),
      other_celltypes_df %>% mutate(group = "Non-PT Cells")
    )
    
    p_comparison <- ggplot(comparison_data, aes(x = cortical_k2f, y = pseudotime, fill = cortical_k2f)) +
      geom_violin(alpha = 0.5) +
      geom_boxplot(width = 0.3, alpha = 0.7) +
      facet_wrap(~group) +
      labs(title = "Cortical K2/F Difference in Pseudotime: PT vs Non-PT Cells",
           x = "Cortical K2/F",
           y = "Pseudotime") +
      theme_classic() +
      theme(legend.position = "none")
    
    print(p_comparison)
    ggsave("fig16_pt_vs_nonpt_comparison_k2f.pdf", p_comparison, width = 10, height = 5)
    ggsave("fig16_pt_vs_nonpt_comparison_k2f.png", p_comparison, width = 10, height = 5, dpi = 300)
  }
}

# 7. Functional enrichment of all PT Cortical K2/F-differential genes
if(tradeseq_pt_success && exists("k2f_results_pt")) {
  sig_pt_genes <- k2f_results_pt %>%
    filter(significant) %>%
    pull(gene)
  
  if(length(sig_pt_genes) > 10) {
    print("\nRunning GO enrichment for all PT genes...")
    
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
        print("GO Enrichment for all PT Cortical K2/F-differential genes:")
        p_go_pt <- dotplot(go_results_pt, showCategory = 20) +
          ggtitle("GO Enrichment: All PT Cortical K2/F-Differential Genes")
        print(p_go_pt)
        ggsave("fig17_go_enrichment_pt_all_k2f.pdf", p_go_pt, width = 10, height = 8)
        ggsave("fig17_go_enrichment_pt_all_k2f.png", p_go_pt, width = 10, height = 8, dpi = 300)
        
        # Save GO results
        write.csv(as.data.frame(go_results_pt), "go_enrichment_pt_all_k2f.csv", row.names = FALSE)
      }
      
      # KEGG pathway enrichment
      kegg_results_pt <- enrichKEGG(
        gene = gene_entrez_pt$ENTREZID,
        organism = "hsa",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05
      )
      
      if(nrow(kegg_results_pt) > 0) {
        print("KEGG Pathway Enrichment for all PT Cortical K2/F-differential genes:")
        p_kegg_pt <- dotplot(kegg_results_pt, showCategory = 15) +
          ggtitle("KEGG Pathways: All PT Cortical K2/F-Differential Genes")
        print(p_kegg_pt)
        ggsave("fig18_kegg_enrichment_pt_all_k2f.pdf", p_kegg_pt, width = 10, height = 8)
        ggsave("fig18_kegg_enrichment_pt_all_k2f.png", p_kegg_pt, width = 10, height = 8, dpi = 300)
        
        # Save KEGG results
        write.csv(as.data.frame(kegg_results_pt), "kegg_enrichment_pt_all_k2f.csv", row.names = FALSE)
      }
    }, error = function(e) {
      print("Error in enrichment for all PT genes:")
      print(e)
    })
  }
}

# ============================================
# PART 7: Functional Enrichment Analysis (Overall)
# ============================================

print("\n========================================")
print("PART 7: FUNCTIONAL ENRICHMENT (OVERALL)")
print("========================================\n")

# Get significant genes from different analyses
if(phenopath_success && nrow(phenopath_results) > 0) {
  sig_phenopath_genes <- phenopath_results %>%
    filter(sig_k2f_effect) %>%
    pull(gene)
} else {
  sig_phenopath_genes <- character(0)
}

sig_tradeseq_genes <- k2f_trajectory_results %>%
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
      p_go_pheno <- dotplot(go_results_pheno, showCategory = 15) +
        ggtitle("GO Enrichment: PhenoPath Cortical K2/F-Interaction Genes")
      print(p_go_pheno)
      ggsave("fig19_go_enrichment_phenopath_k2f.pdf", p_go_pheno, width = 10, height = 8)
      ggsave("fig19_go_enrichment_phenopath_k2f.png", p_go_pheno, width = 10, height = 8, dpi = 300)
      
      write.csv(as.data.frame(go_results_pheno), "go_enrichment_phenopath_k2f.csv", row.names = FALSE)
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
      p_go_trade <- dotplot(go_results_trade, showCategory = 15) +
        ggtitle("GO Enrichment: TradeSeq Cortical K2/F-Differential Genes")
      print(p_go_trade)
      ggsave("fig20_go_enrichment_tradeseq_k2f.pdf", p_go_trade, width = 10, height = 8)
      ggsave("fig20_go_enrichment_tradeseq_k2f.png", p_go_trade, width = 10, height = 8, dpi = 300)
      
      write.csv(as.data.frame(go_results_trade), "go_enrichment_tradeseq_k2f.csv", row.names = FALSE)
    }
  }, error = function(e) {
    print("Error in GO enrichment for tradeSeq genes:")
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
  cat(paste("   - Genes with strong Cortical K2/F effects:", sum(phenopath_results$sig_k2f_effect), "\n"))
  if(exists("cor_overall")) {
    cat(paste("   - Correlation with Slingshot:", round(cor_overall, 3), "\n\n"))
  }
} else {
  cat("2. PHENOPATH RESULTS:\n")
  cat("   - PhenoPath analysis failed or was skipped\n\n")
}

cat("3. TRADESEQ RESULTS (ALL CELLS):\n")
cat(paste("   - Significant genes (padj < 0.05):", sum(k2f_trajectory_results$significant), "\n"))
cat(paste("   - Significant genes (padj < 0.01):", sum(k2f_trajectory_results$padj < 0.01), "\n\n"))

if(phenopath_success && exists("combined_results")) {
  cat("4. METHOD COMPARISON:\n")
  cat(paste("   - Genes identified by PhenoPath:", sum(phenopath_results$sig_k2f_effect), "\n"))
  cat(paste("   - Genes identified by tradeSeq:", sum(k2f_trajectory_results$significant), "\n"))
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

cat("5. ALL PT CELLS SPECIFIC ANALYSIS:\n")
cat(paste("   - All PT Cells Wilcoxon p-value:", signif(pt_all_test$p.value, 3), "\n"))
cat(paste("   - All PT Cells Cohen's d:", round(pt_all_cohens$estimate, 3), "\n"))
if(tradeseq_pt_success) {
  cat(paste("   - All PT Cortical K2/F-differential genes:", sum(k2f_results_pt$significant), "\n\n"))
} else {
  cat("   - TradeSeq analysis skipped for PT cells\n\n")
}

cat("OUTPUT FILES CREATED:\n")
cat("   CSV FILES:\n")
cat("   - trajectory_stats_cortical_k2f.csv\n")
cat("   - celltype_composition_by_k2f.csv\n")
if(phenopath_success) {
  cat("   - phenopath_cortical_k2f_interaction_genes.csv\n")
  cat("   - phenopath_slingshot_comparison.csv\n")
  cat("   - correlation_by_k2f_level.csv\n")
}
cat("   - tradeseq_cortical_k2f_specific_trajectory_genes.csv\n")
if(phenopath_success && exists("combined_results")) {
  cat("   - combined_phenopath_tradeseq_cortical_k2f_results.csv\n")
}
cat("   - pt_all_cells_stats_k2f.csv\n")
if(tradeseq_pt_success) {
  cat("   - pt_all_cells_cortical_k2f_specific_genes.csv\n")
  cat("   - pt_all_expr_summary_k2f.csv\n")
}
if(exists("go_results_pt") && nrow(go_results_pt) > 0) {
  cat("   - go_enrichment_pt_all_k2f.csv\n")
}
if(exists("kegg_results_pt") && nrow(kegg_results_pt) > 0) {
  cat("   - kegg_enrichment_pt_all_k2f.csv\n")
}
if(exists("go_results_pheno") && nrow(go_results_pheno) > 0) {
  cat("   - go_enrichment_phenopath_k2f.csv\n")
}
if(exists("go_results_trade") && nrow(go_results_trade) > 0) {
  cat("   - go_enrichment_tradeseq_k2f.csv\n")
}

cat("\n   FIGURE FILES (PDF and PNG):\n")
cat("   - fig1_trajectory_comparison_cortical_k2f\n")
cat("   - fig2_celltype_composition_cortical_k2f\n")
if(phenopath_success) {
  cat("   - fig3_slingshot_vs_phenopath_comparison\n")
  cat("   - fig4_pseudotime_umap_comparison\n")
  cat("   - fig5_top_phenopath_genes\n")
  cat("   - fig6_phenopath_vs_tradeseq_comparison\n")
  cat("   - fig7_venn_diagram_methods_overlap\n")
  cat("   - fig8_top_genes_both_methods\n")
}
cat("   - fig9_pt_all_density_k2f\n")
cat("   - fig10_pt_all_boxplot_k2f\n")
cat("   - fig11_pt_all_umap_k2f\n")
cat("   - fig12_pt_all_ridge_k2f\n")
cat("   - fig13_pt_by_celltype_k2f\n")
if(tradeseq_pt_success) {
  cat("   - fig14_top_pt_all_k2f_genes\n")
  cat("   - fig15_pt_all_expr_pattern_k2f\n")
}
if(exists("p_comparison")) {
  cat("   - fig16_pt_vs_nonpt_comparison_k2f\n")
}
if(exists("go_results_pt") && nrow(go_results_pt) > 0) {
  cat("   - fig17_go_enrichment_pt_all_k2f\n")
}
if(exists("kegg_results_pt") && nrow(kegg_results_pt) > 0) {
  cat("   - fig18_kegg_enrichment_pt_all_k2f\n")
}
if(exists("go_results_pheno") && nrow(go_results_pheno) > 0) {
  cat("   - fig19_go_enrichment_phenopath_k2f\n")
}
if(exists("go_results_trade") && nrow(go_results_trade) > 0) {
  cat("   - fig20_go_enrichment_tradeseq_k2f\n")
}

cat("\n========================================\n")
cat("Analysis complete!\n")
cat("========================================\n\n")

# ============================================
# ADDITIONAL: Create summary statistics table
# ============================================

summary_stats <- data.frame(
  Analysis = c("Overall Trajectory", 
               "PhenoPath", 
               "TradeSeq (All Cells)", 
               "All PT Cells"),
  Wilcoxon_pvalue = c(
    signif(wilcox_result$p.value, 3),
    NA,
    NA,
    signif(pt_all_test$p.value, 3)
  ),
  Cohens_d = c(
    round(cohens_d$estimate, 3),
    NA,
    NA,
    round(pt_all_cohens$estimate, 3)
  ),
  Significant_genes = c(
    NA,
    ifelse(phenopath_success, sum(phenopath_results$sig_k2f_effect), NA),
    sum(k2f_trajectory_results$significant),
    ifelse(tradeseq_pt_success, sum(k2f_results_pt$significant), NA)
  )
)

write.csv(summary_stats, "summary_statistics_cortical_k2f.csv", row.names = FALSE)
print("\nSummary Statistics Table:")
print(summary_stats)

# ============================================
# SAVE SESSION INFO
# ============================================

session_info <- sessionInfo()
sink("session_info.txt")
print(session_info)
sink()

cat("\nSession info saved to session_info.txt\n")
cat("\n========================================\n")
cat("ALL ANALYSES AND PLOTS COMPLETE!\n")
cat("========================================\n")









