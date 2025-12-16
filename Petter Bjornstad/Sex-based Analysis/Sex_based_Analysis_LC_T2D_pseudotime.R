#Sex based analysis pseudotime LC & T2D 

# Load required libraries (already loaded in your scripts)
library(Seurat)
library(dplyr)
library(tidyverse)

# ============================================================================
# STEP 1: Load T2D Seurat Object
# ============================================================================
load('C:/Users/netio/Documents/UofW/Rockies/ROCKIES_T2D_SGLT2_DylanEdits_Line728.RData')
so_t2d <- so_kpmp_sc
remove(so_kpmp_sc)

# Apply celltype classifications
so_t2d$celltype1 <- case_when(grepl("PT-",so_t2d$celltype_rpca)~"PT",
                              grepl("TAL-",so_t2d$celltype_rpca)~"TAL",
                              grepl("EC-",so_t2d$celltype_rpca)~"EC",
                              grepl("POD",so_t2d$celltype_rpca)~"POD",
                              grepl("MAC",so_t2d$celltype_rpca)~"MAC",
                              grepl("MON",so_t2d$celltype_rpca)~"MON",
                              grepl("PC-",so_t2d$celltype_rpca)~"PC",
                              grepl("FIB",so_t2d$celltype_rpca)~"FIB_MC_VSMC",
                              grepl("DTL",so_t2d$celltype_rpca)~"DTL",
                              so_t2d$celltype_rpca=="DCT"~"DCT",
                              so_t2d$celltype_rpca=="ATL"~"ATL",
                              so_t2d$celltype_rpca=="B"~"B",
                              so_t2d$celltype_rpca=="T"~"T")
so_t2d$celltype1 <- as.character(so_t2d$celltype1)
so_t2d$KPMP_celltype2 <- as.character(so_t2d$KPMP_celltype)
so_t2d$celltype2 <- ifelse(so_t2d$KPMP_celltype=="aPT" | 
                             so_t2d$KPMP_celltype=="PT-S1/S2" | 
                             so_t2d$KPMP_celltype == "PT-S3","PT",
                           ifelse(grepl("TAL",so_t2d$KPMP_celltype),"TAL",
                                  ifelse(grepl("EC-",so_t2d$KPMP_celltype),"EC",so_t2d$KPMP_celltype2)))
so_t2d$DCT_celltype <- ifelse((so_t2d$KPMP_celltype=="DCT" | 
                                 so_t2d$KPMP_celltype=="dDCT"), "DCT","Non-DCT")

# Filter T2D samples
so_t2d <- subset(so_t2d, subset = record_id != 'CRC-55')
so_t2d <- subset(so_t2d, subset = group == 'Type_2_Diabetes')

cat("T2D object:\n")
cat("  Cells:", ncol(so_t2d), "\n")
cat("  Samples:", length(unique(so_t2d$record_id)), "\n")
print(table(so_t2d$record_id))

# ============================================================================
# STEP 2: Load Lean Control Seurat Object
# ============================================================================
load('C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/No_Med_line700.Rdata')
so_lean <- so_kpmp_sc
remove(so_kpmp_sc)

# Apply same celltype classifications to Lean Control object
so_lean$celltype1 <- case_when(grepl("PT-",so_lean$celltype_rpca)~"PT",
                               grepl("TAL-",so_lean$celltype_rpca)~"TAL",
                               grepl("EC-",so_lean$celltype_rpca)~"EC",
                               grepl("POD",so_lean$celltype_rpca)~"POD",
                               grepl("MAC",so_lean$celltype_rpca)~"MAC",
                               grepl("MON",so_lean$celltype_rpca)~"MON",
                               grepl("PC-",so_lean$celltype_rpca)~"PC",
                               grepl("FIB",so_lean$celltype_rpca)~"FIB_MC_VSMC",
                               grepl("DTL",so_lean$celltype_rpca)~"DTL",
                               so_lean$celltype_rpca=="DCT"~"DCT",
                               so_lean$celltype_rpca=="ATL"~"ATL",
                               so_lean$celltype_rpca=="B"~"B",
                               so_lean$celltype_rpca=="T"~"T")
so_lean$celltype1 <- as.character(so_lean$celltype1)
so_lean$KPMP_celltype2 <- as.character(so_lean$KPMP_celltype)
so_lean$celltype2 <- ifelse(so_lean$KPMP_celltype=="aPT" | 
                              so_lean$KPMP_celltype=="PT-S1/S2" | 
                              so_lean$KPMP_celltype == "PT-S3","PT",
                            ifelse(grepl("TAL",so_lean$KPMP_celltype),"TAL",
                                   ifelse(grepl("EC-",so_lean$KPMP_celltype),"EC",so_lean$KPMP_celltype2)))
so_lean$DCT_celltype <- ifelse((so_lean$KPMP_celltype=="DCT" | 
                                  so_lean$KPMP_celltype=="dDCT"), "DCT","Non-DCT")

# Filter Lean Control samples
so_lean <- subset(so_lean, subset = record_id != 'CRC-55')
so_lean <- subset(so_lean, subset = group == 'Lean_Control')

cat("\nLean Control object:\n")
cat("  Cells:", ncol(so_lean), "\n")
cat("  Samples:", length(unique(so_lean$record_id)), "\n")
print(table(so_lean$record_id))

# ============================================================================
# STEP 3: Check for Overlapping Samples
# ============================================================================
t2d_samples <- unique(so_t2d$record_id)
lean_samples <- unique(so_lean$record_id)
overlapping_samples <- intersect(t2d_samples, lean_samples)

cat("\nOverlapping samples between T2D and Lean Control:", length(overlapping_samples), "\n")
if(length(overlapping_samples) > 0) {
  cat("Overlapping sample IDs:\n")
  print(overlapping_samples)
}

# ============================================================================
# STEP 4: Merge Seurat Objects (removing duplicate samples if any)
# ============================================================================
if(length(overlapping_samples) > 0) {
  # Remove overlapping samples from Lean Control object
  cat("\nRemoving", length(overlapping_samples), "overlapping samples from Lean Control object...\n")
  so_lean <- subset(so_lean, subset = !record_id %in% overlapping_samples)
}

# Merge the objects
so_merged <- merge(x = so_t2d, 
                   y = so_lean,
                   add.cell.ids = c("T2D", "Lean"),
                   project = "T2D_vs_Lean_Sex_Analysis")

cat("\nMerged object:\n")
cat("  Total cells:", ncol(so_merged), "\n")
cat("  Total samples:", length(unique(so_merged$record_id)), "\n")
cat("  T2D samples:", sum(so_merged$group == "Type_2_Diabetes" & 
                            !duplicated(so_merged$record_id)), "\n")
cat("  Lean Control samples:", sum(so_merged$group == "Lean_Control" & 
                                     !duplicated(so_merged$record_id)), "\n")

# ============================================================================
# STEP 5: Verify merged object and check for sex information
# ============================================================================
# Check sex distribution
cat("\nSex distribution:\n")
if("sex" %in% colnames(so_merged@meta.data)) {
  print(table(so_merged$group, so_merged$sex))
} else if("Sex" %in% colnames(so_merged@meta.data)) {
  print(table(so_merged$group, so_merged$Sex))
} else {
  cat("No 'sex' or 'Sex' column found in metadata. Available columns:\n")
  print(colnames(so_merged@meta.data))
}

# Summary by group and sample
summary_df <- so_merged@meta.data %>%
  dplyr::select(record_id, group) %>%
  filter(!duplicated(record_id)) %>%
  group_by(group) %>%
  summarise(n_samples = n())

cat("\nSummary by group:\n")
print(summary_df)

# ============================================================================
# STEP 6: Save merged object
# ============================================================================
dir.results <- 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/'
dir.create(dir.results, recursive = TRUE, showWarnings = FALSE)

save(so_merged, file = paste0(dir.results, 'T2D_Lean_Merged_Seurat.RData'))
cat("\nMerged Seurat object saved to:", paste0(dir.results, 'T2D_Lean_Merged_Seurat.RData'), "\n")

# Clean up individual objects to free memory
remove(so_t2d, so_lean)
gc()

cat("\nMerge complete! Use 'so_merged' for your sex and T2D status analysis.\n")









###Load data














##Analyze Data 

########## Lineage Tracing (Slingshot) Analysis - T2D vs Lean Control (PT Cells)
########## Colored by Sex and Diabetes Status

library(slingshot)
library(ggplot2)
library(viridis)
library(ggridges)
library(patchwork)
library(dplyr)
library(ggpubr)
library(rstatix)
library(Seurat)

# ============================================================================
# Load merged data (assuming you've already created so_merged from previous step)
# ============================================================================
# If not already loaded:
# load('C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Lean_Merged_Seurat.RData')

dir.results <- 'C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis/T2D_Lean_Merged/pseudotime/'
dir.create(dir.results, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Subset to PT cells and run PCA + UMAP
# ============================================================================
so_pt <- subset(so_merged, subset = celltype2 == 'PT')

cat("PT cells selected:\n")
cat("  Total cells:", ncol(so_pt), "\n")
cat("  By group:\n")
print(table(so_pt$group))
cat("  By sex:\n")
print(table(so_pt$sex))
cat("  By group and sex:\n")
print(table(so_pt$group, so_pt$sex))

# Join layers if using Seurat v5 (before normalization)
if ("layers" %in% slotNames(so_pt[["RNA"]])) {
  cat("\nJoining Seurat v5 layers...\n")
  so_pt[["RNA"]] <- JoinLayers(so_pt[["RNA"]])
}

# Normalize and scale data
cat("Normalizing data...\n")
so_pt <- NormalizeData(so_pt)
so_pt <- FindVariableFeatures(so_pt, selection.method = "vst", nfeatures = 2000)
so_pt <- ScaleData(so_pt)

# Run PCA
cat("Running PCA...\n")
so_pt <- RunPCA(so_pt, features = VariableFeatures(object = so_pt))

# Check PCA - optional
ElbowPlot(so_pt, ndims = 50)

# Run UMAP
cat("Running UMAP...\n")
so_pt <- RunUMAP(so_pt, dims = 1:30)

cat("\nPCA and UMAP complete.\n")

# ============================================================================
# Run Slingshot
# ============================================================================
cat("Running Slingshot pseudotime analysis...\n")
sling_res <- slingshot(as.SingleCellExperiment(so_pt), 
                       clusterLabels = 'KPMP_celltype', 
                       start.clus = 'PT-S1', 
                       end.clus = 'aPT', 
                       reducedDim = 'UMAP')

so_pt$pseudotime <- slingPseudotime(sling_res)[,1]

cat("\nPseudotime calculated. NA values:", sum(is.na(so_pt$pseudotime)), "\n")

# ============================================================================
# Prepare plotting data
# ============================================================================
umap_coords <- Embeddings(so_pt, reduction = "umap")
sling_curves <- slingCurves(sling_res)

plot_df <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  pseudotime = so_pt$pseudotime,
  celltype = so_pt$KPMP_celltype,
  sex = so_pt$sex,
  group = so_pt$group
)

# Clean data (remove NA pseudotime)
plot_df_clean <- plot_df %>% filter(!is.na(pseudotime))

# Create cleaner labels for plotting
plot_df_clean <- plot_df_clean %>%
  mutate(
    group_label = case_when(
      group == "Type_2_Diabetes" ~ "T2D",
      group == "Lean_Control" ~ "Lean Control",
      TRUE ~ as.character(group)
    ),
    sex_group = paste(sex, group_label, sep = " - ")
  )

cat("\nPlotting data prepared. Creating visualizations...\n")

# ============================================================================
# Plot 1: UMAP with pseudotime
# ============================================================================
cat("Creating Plot 1: UMAP with pseudotime...\n")
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

# ============================================================================
# Plot 2: UMAP colored by Sex
# ============================================================================
cat("Creating Plot 2: UMAP colored by Sex...\n")
p2_sex <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = sex)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = c("Female" = "#E64B35", "Male" = "#4DBBD5")) +
  theme_classic() +
  labs(title = "UMAP Colored by Sex",
       color = "Sex") +
  theme(legend.position = "right")

for(i in seq_along(sling_curves)) {
  curve_coords <- sling_curves[[i]]$s[sling_curves[[i]]$ord, ]
  p2_sex <- p2_sex + geom_path(data = data.frame(UMAP_1 = curve_coords[, 1], 
                                                 UMAP_2 = curve_coords[, 2]),
                               aes(x = UMAP_1, y = UMAP_2),
                               color = "black", size = 2, inherit.aes = FALSE)
}

print(p2_sex)
ggsave(paste0(dir.results, "Slingshot_UMAP_by_Sex.pdf"), 
       plot = p2_sex, width = 8, height = 6)

# ============================================================================
# Plot 3: UMAP colored by Diabetes Status (Group)
# ============================================================================
cat("Creating Plot 3: UMAP colored by Diabetes Status...\n")
p3_group <- ggplot(plot_df_clean, aes(x = UMAP_1, y = UMAP_2, color = group_label)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = c("T2D" = "#DC0000", "Lean Control" = "#00A087")) +
  theme_classic() +
  labs(title = "UMAP Colored by Diabetes Status",
       color = "Group") +
  theme(legend.position = "right")

for(i in seq_along(sling_curves)) {
  curve_coords <- sling_curves[[i]]$s[sling_curves[[i]]$ord, ]
  p3_group <- p3_group + geom_path(data = data.frame(UMAP_1 = curve_coords[, 1], 
                                                     UMAP_2 = curve_coords[, 2]),
                                   aes(x = UMAP_1, y = UMAP_2),
                                   color = "black", size = 2, inherit.aes = FALSE)
}

print(p3_group)
ggsave(paste0(dir.results, "Slingshot_UMAP_by_Group.pdf"), 
       plot = p3_group, width = 8, height = 6)

# ============================================================================
# Plot 4: UMAP split by Sex (faceted)
# ============================================================================
cat("Creating Plot 4: UMAP faceted by Sex...\n")
p4_sex_facet <- ggplot(plot_df_clean, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_viridis(option = "plasma") +
  facet_wrap(~sex) +
  theme_classic() +
  labs(title = "Slingshot Pseudotime by Sex",
       color = "Pseudotime")

for(i in seq_along(sling_curves)) {
  curve_coords <- sling_curves[[i]]$s[sling_curves[[i]]$ord, ]
  p4_sex_facet <- p4_sex_facet + geom_path(data = data.frame(UMAP_1 = curve_coords[, 1], 
                                                             UMAP_2 = curve_coords[, 2]),
                                           aes(x = UMAP_1, y = UMAP_2),
                                           color = "black", size = 1.5, inherit.aes = FALSE)
}

print(p4_sex_facet)
ggsave(paste0(dir.results, "Slingshot_UMAP_Facet_by_Sex.pdf"), 
       plot = p4_sex_facet, width = 12, height = 5)

# ============================================================================
# Plot 5: UMAP split by Diabetes Status (faceted)
# ============================================================================
cat("Creating Plot 5: UMAP faceted by Diabetes Status...\n")
p5_group_facet <- ggplot(plot_df_clean, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_viridis(option = "plasma") +
  facet_wrap(~group_label) +
  theme_classic() +
  labs(title = "Slingshot Pseudotime by Diabetes Status",
       color = "Pseudotime")

for(i in seq_along(sling_curves)) {
  curve_coords <- sling_curves[[i]]$s[sling_curves[[i]]$ord, ]
  p5_group_facet <- p5_group_facet + geom_path(data = data.frame(UMAP_1 = curve_coords[, 1], 
                                                                 UMAP_2 = curve_coords[, 2]),
                                               aes(x = UMAP_1, y = UMAP_2),
                                               color = "black", size = 1.5, inherit.aes = FALSE)
}

print(p5_group_facet)
ggsave(paste0(dir.results, "Slingshot_UMAP_Facet_by_Group.pdf"), 
       plot = p5_group_facet, width = 12, height = 5)

# ============================================================================
# Plot 6: UMAP split by both Sex AND Diabetes Status (2x2 facet)
# ============================================================================
cat("Creating Plot 6: UMAP faceted by Sex and Diabetes Status...\n")
p6_both_facet <- ggplot(plot_df_clean, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_viridis(option = "plasma") +
  facet_grid(sex ~ group_label) +
  theme_classic() +
  labs(title = "Slingshot Pseudotime by Sex and Diabetes Status",
       color = "Pseudotime")

for(i in seq_along(sling_curves)) {
  curve_coords <- sling_curves[[i]]$s[sling_curves[[i]]$ord, ]
  p6_both_facet <- p6_both_facet + geom_path(data = data.frame(UMAP_1 = curve_coords[, 1], 
                                                               UMAP_2 = curve_coords[, 2]),
                                             aes(x = UMAP_1, y = UMAP_2),
                                             color = "black", size = 1.5, inherit.aes = FALSE)
}

print(p6_both_facet)
ggsave(paste0(dir.results, "Slingshot_UMAP_Facet_Sex_and_Group.pdf"), 
       plot = p6_both_facet, width = 12, height = 10)

# ============================================================================
# Summary Statistics
# ============================================================================
cat("Calculating summary statistics...\n")
pseudotime_summary <- plot_df_clean %>%
  group_by(sex, group_label) %>%
  summarise(
    n = n(),
    mean_pseudotime = mean(pseudotime),
    median_pseudotime = median(pseudotime),
    sd_pseudotime = sd(pseudotime),
    se_pseudotime = sd(pseudotime) / sqrt(n()),
    .groups = 'drop'
  )

print(pseudotime_summary)
write.csv(pseudotime_summary, paste0(dir.results, "Pseudotime_Summary_by_Sex_and_Group.csv"), 
          row.names = FALSE)

# ============================================================================
# Plot 7: Violin plot - Pseudotime by Sex
# ============================================================================
cat("Creating Plot 7: Violin plot by Sex...\n")
stat_test_sex <- plot_df_clean %>%
  wilcox_test(pseudotime ~ sex) %>%
  add_significance()

print("Statistical test - Sex:")
print(stat_test_sex)

p7_violin_sex <- ggplot(plot_df_clean, aes(x = sex, y = pseudotime, fill = sex)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = c("Female" = "#E64B35", "Male" = "#4DBBD5")) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_classic() +
  labs(title = "Pseudotime Comparison by Sex",
       x = "Sex",
       y = "Pseudotime") +
  theme(legend.position = "none")

print(p7_violin_sex)
ggsave(paste0(dir.results, "Pseudotime_by_Sex_Violin.pdf"), 
       plot = p7_violin_sex, width = 6, height = 6)

# ============================================================================
# Plot 8: Violin plot - Pseudotime by Diabetes Status
# ============================================================================
cat("Creating Plot 8: Violin plot by Diabetes Status...\n")
stat_test_group <- plot_df_clean %>%
  wilcox_test(pseudotime ~ group_label) %>%
  add_significance()

print("Statistical test - Diabetes Status:")
print(stat_test_group)

p8_violin_group <- ggplot(plot_df_clean, aes(x = group_label, y = pseudotime, fill = group_label)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = c("T2D" = "#DC0000", "Lean Control" = "#00A087")) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_classic() +
  labs(title = "Pseudotime Comparison by Diabetes Status",
       x = "Group",
       y = "Pseudotime") +
  theme(legend.position = "none")

print(p8_violin_group)
ggsave(paste0(dir.results, "Pseudotime_by_Group_Violin.pdf"), 
       plot = p8_violin_group, width = 6, height = 6)

# ============================================================================
# Plot 9: Violin plot - Pseudotime by Sex AND Diabetes Status
# ============================================================================
cat("Creating Plot 9: Violin plot by Sex and Diabetes Status...\n")
# Statistical tests for each group
stat_test_combined <- plot_df_clean %>%
  group_by(group_label) %>%
  wilcox_test(pseudotime ~ sex) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "sex", dodge = 0.8)

print("Statistical test - Sex within each Diabetes Status:")
print(stat_test_combined)

p9_violin_combined <- ggplot(plot_df_clean, aes(x = sex, y = pseudotime, fill = group_label)) +
  geom_violin(alpha = 0.6, trim = FALSE, position = position_dodge(0.8)) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("T2D" = "#DC0000", "Lean Control" = "#00A087")) +
  stat_pvalue_manual(stat_test_combined, label = "p.adj.signif", 
                     tip.length = 0.01, hide.ns = FALSE) +
  theme_classic() +
  labs(title = "Pseudotime Comparison by Sex and Diabetes Status",
       x = "Sex",
       y = "Pseudotime",
       fill = "Group") +
  theme(legend.position = "right")

print(p9_violin_combined)
ggsave(paste0(dir.results, "Pseudotime_by_Sex_and_Group_Violin.pdf"), 
       plot = p9_violin_combined, width = 8, height = 6)

# ============================================================================
# Plot 10: Density plot - Pseudotime by Sex
# ============================================================================
cat("Creating Plot 10: Density plot by Sex...\n")
p10_density_sex <- ggplot(plot_df_clean, aes(x = pseudotime, fill = sex)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Female" = "#E64B35", "Male" = "#4DBBD5")) +
  theme_classic() +
  labs(title = "Pseudotime Distribution by Sex",
       x = "Pseudotime",
       y = "Density",
       fill = "Sex")

print(p10_density_sex)
ggsave(paste0(dir.results, "Pseudotime_Density_by_Sex.pdf"), 
       plot = p10_density_sex, width = 8, height = 6)

# ============================================================================
# Plot 11: Density plot - Pseudotime by Diabetes Status
# ============================================================================
cat("Creating Plot 11: Density plot by Diabetes Status...\n")
p11_density_group <- ggplot(plot_df_clean, aes(x = pseudotime, fill = group_label)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("T2D" = "#DC0000", "Lean Control" = "#00A087")) +
  theme_classic() +
  labs(title = "Pseudotime Distribution by Diabetes Status",
       x = "Pseudotime",
       y = "Density",
       fill = "Group")

print(p11_density_group)
ggsave(paste0(dir.results, "Pseudotime_Density_by_Group.pdf"), 
       plot = p11_density_group, width = 8, height = 6)

# ============================================================================
# Plot 12: Density plot - Pseudotime by Sex within each Diabetes Status
# ============================================================================
cat("Creating Plot 12: Density plot by Sex within Diabetes Status...\n")
p12_density_combined <- ggplot(plot_df_clean, aes(x = pseudotime, fill = sex)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Female" = "#E64B35", "Male" = "#4DBBD5")) +
  facet_wrap(~group_label) +
  theme_classic() +
  labs(title = "Pseudotime Distribution by Sex within Diabetes Status",
       x = "Pseudotime",
       y = "Density",
       fill = "Sex")

print(p12_density_combined)
ggsave(paste0(dir.results, "Pseudotime_Density_by_Sex_within_Group.pdf"), 
       plot = p12_density_combined, width = 12, height = 5)

# ============================================================================
# Plot 13: Ridgeline plot - Cell type distribution by Sex
# ============================================================================
cat("Creating Plot 13: Ridgeline plot by Sex...\n")
p13_ridge_sex <- ggplot(plot_df_clean, aes(x = pseudotime, y = celltype, fill = sex)) +
  geom_density_ridges(alpha = 0.6, scale = 2) +
  scale_fill_manual(values = c("Female" = "#E64B35", "Male" = "#4DBBD5")) +
  theme_classic() +
  labs(title = "Cell Type Distribution Along Pseudotime by Sex",
       x = "Pseudotime",
       y = "Cell Type",
       fill = "Sex") +
  theme(plot.title = element_text(hjust = 0.5))

print(p13_ridge_sex)
ggsave(paste0(dir.results, "Pseudotime_Ridge_by_Celltype_Sex.pdf"), 
       plot = p13_ridge_sex, width = 10, height = 6)

# ============================================================================
# Plot 14: Ridgeline plot - Cell type distribution by Diabetes Status
# ============================================================================
cat("Creating Plot 14: Ridgeline plot by Diabetes Status...\n")
p14_ridge_group <- ggplot(plot_df_clean, aes(x = pseudotime, y = celltype, fill = group_label)) +
  geom_density_ridges(alpha = 0.6, scale = 2) +
  scale_fill_manual(values = c("T2D" = "#DC0000", "Lean Control" = "#00A087")) +
  theme_classic() +
  labs(title = "Cell Type Distribution Along Pseudotime by Diabetes Status",
       x = "Pseudotime",
       y = "Cell Type",
       fill = "Group") +
  theme(plot.title = element_text(hjust = 0.5))

print(p14_ridge_group)
ggsave(paste0(dir.results, "Pseudotime_Ridge_by_Celltype_Group.pdf"), 
       plot = p14_ridge_group, width = 10, height = 6)

# ============================================================================
# Combined Panel Plots
# ============================================================================
cat("Creating combined panel plots...\n")

# Panel 1: Overview
combined_plot1 <- (p1 | p2_sex | p3_group) / (p7_violin_sex | p8_violin_group | p9_violin_combined)
print(combined_plot1)
ggsave(paste0(dir.results, "Combined_Panel_1_Overview.pdf"), 
       plot = combined_plot1, width = 18, height = 10)

# Panel 2: Faceted UMAPs
combined_plot2 <- (p4_sex_facet / p5_group_facet)
print(combined_plot2)
ggsave(paste0(dir.results, "Combined_Panel_2_Faceted_UMAPs.pdf"), 
       plot = combined_plot2, width = 12, height = 10)

# Panel 3: Density analyses
combined_plot3 <- (p10_density_sex | p11_density_group) / p12_density_combined
print(combined_plot3)
ggsave(paste0(dir.results, "Combined_Panel_3_Densities.pdf"), 
       plot = combined_plot3, width = 14, height = 10)

# Panel 4: Main figure with 2x2 facet
combined_plot4 <- p6_both_facet / (p9_violin_combined | p12_density_combined)
print(combined_plot4)
ggsave(paste0(dir.results, "Combined_Panel_4_Main_Figure.pdf"), 
       plot = combined_plot4, width = 14, height = 12)

# ============================================================================
# Save processed Seurat object with pseudotime
# ============================================================================
cat("Saving processed Seurat object...\n")
save(so_pt, sling_res, file = paste0(dir.results, "PT_Cells_with_Pseudotime.RData"))

cat("\n=== Analysis Complete ===\n")
cat("Results saved to:", dir.results, "\n")
cat("\nKey findings:\n")
print(pseudotime_summary)
cat("\nStatistical tests saved to CSV files in the results directory.\n")
















