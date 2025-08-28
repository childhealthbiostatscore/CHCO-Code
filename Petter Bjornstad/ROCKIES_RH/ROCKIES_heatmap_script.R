###Heatmap Analysis 

#packages 






















#heatmaps 





total_results <- data.table::fread(paste0(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_TCA_cycle_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
heatmap_data$Variable <- factor(heatmap_data$Variable, levels = custom_order)
custom_labels <- c("Average Cortical K2","Average Medulla K2","Average Cortical F","Average Medulla F",
                   "Average Cortical K2/F","Average Medulla K2/F")

heat_map_p <- ggplot(heatmap_data, aes(x = Variable, y = Gene, fill = logFC)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = signif), size = 3, color = "black") +
  scale_fill_gradient2(
    low = "#264653", mid = "white", high = "darkred",
    midpoint = 0,
    name = "LogFC"
  ) +
  theme_minimal() +
  labs(title = "TCA Cycle Genes vs. PET Variables (T2D)",
       subtitle = "Proximal Tubule Cells",
       x = "Exposure",
       y = "Gene") +
  scale_x_discrete(labels = setNames(custom_labels, custom_order))+
  theme(
    text = element_text(face="bold"),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  )

# custom_colors <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51", "darkred")

png(fs::path(dir.results, "Heatmap_TCA_cycle_NEBULA_PT_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()






total_results <- data.table::fread(paste0(dir.results,"NEBULA_ox_phos_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))
# total_results <- read.csv(fs::path(dir.results,"NEBULA_ox_phos_PT_cells_PET_Variables_unadjusted_pooled_offset.csv"))

# Define significance stars
total_results <- total_results %>%
  mutate(signif = case_when(
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Select only the needed columns and rename LogFC for clarity
heatmap_data <- total_results %>%
  dplyr::select(Gene, Variable, logFC, signif)

custom_order <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
heatmap_data$Variable <- factor(heatmap_data$Variable, levels = custom_order)
custom_labels <- c("Average Cortical K2","Average Medulla K2","Average Cortical F","Average Medulla F",
                   "Average Cortical K2/F","Average Medulla K2/F")

heat_map_p <- ggplot(heatmap_data, aes(x = Variable, y = Gene, fill = logFC)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = signif), size = 3, color = "black") +
  scale_fill_gradient2(
    low = "#264653", mid = "white", high = "darkred",
    midpoint = 0,
    name = "LogFC"
  ) +
  theme_minimal() +
  labs(title = "Ox Phos Genes vs. PET Variables (T2D)",
       subtitle = "Proximal Tubule Cells",
       x = "Exposure",
       y = "Gene") +
  scale_x_discrete(labels = setNames(custom_labels, custom_order))+
  theme(
    text = element_text(face="bold"),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  )

# custom_colors <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51", "darkred")

png(fs::path(dir.results, "Heatmap_ox_phos_NEBULA_PT_PET_unadjusted_pooled_offset_T2D.png"), 
    width = 1500, height = 2000, res = 300)
print(heat_map_p)
dev.off()


































