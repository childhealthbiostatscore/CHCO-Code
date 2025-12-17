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
                       start.clus = 'PT-S1/S2', 
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











##Interaction analysis 



# ============================================================================
# Statistical Analysis: Testing Group Differences and Interactions
# ============================================================================
cat("\n=== Statistical Analysis of Pseudotime ===\n")

# Prepare data for modeling
model_data <- plot_df_clean %>%
  mutate(
    sex = factor(sex),
    group_label = factor(group_label, levels = c("Lean Control", "T2D"))
  )

# ============================================================================
# 1. Linear Model: Main Effects and Interaction
# ============================================================================
cat("\n--- Linear Model Analysis ---\n")

# Fit linear model with interaction
lm_model <- lm(pseudotime ~ sex * group_label, data = model_data)

# Summary
cat("\nModel Summary:\n")
print(summary(lm_model))

# ANOVA to test significance of terms
cat("\nANOVA Table:\n")
anova_results <- anova(lm_model)
print(anova_results)

# Save results
sink(paste0(dir.results, "Linear_Model_Results.txt"))
cat("Linear Model: pseudotime ~ sex * group_label\n")
cat("==============================================\n\n")
print(summary(lm_model))
cat("\n\nANOVA Table:\n")
print(anova_results)
sink()

# Extract coefficients
coef_summary <- summary(lm_model)$coefficients
write.csv(coef_summary, paste0(dir.results, "Linear_Model_Coefficients.csv"))

# ============================================================================
# 2. Effect Size Calculations
# ============================================================================
cat("\n--- Effect Sizes ---\n")

# Calculate marginal means
library(emmeans)
emm <- emmeans(lm_model, ~ sex * group_label)
cat("\nEstimated Marginal Means:\n")
print(emm)

# Pairwise comparisons
pairs_all <- pairs(emm, adjust = "bonferroni")
cat("\nPairwise Comparisons (Bonferroni corrected):\n")
print(pairs_all)

# Simple effects: Test sex effect within each group
cat("\nSimple Effects - Sex within each Diabetes Status:\n")
simple_effects_group <- emm %>%
  contrast("pairwise", by = "group_label", adjust = "bonferroni")
print(simple_effects_group)

# Simple effects: Test group effect within each sex
cat("\nSimple Effects - Diabetes Status within each Sex:\n")
simple_effects_sex <- emm %>%
  contrast("pairwise", by = "sex", adjust = "bonferroni")
print(simple_effects_sex)

# Save emmeans results
write.csv(as.data.frame(emm), paste0(dir.results, "Estimated_Marginal_Means.csv"))
write.csv(as.data.frame(pairs_all), paste0(dir.results, "Pairwise_Comparisons_All.csv"))
write.csv(as.data.frame(simple_effects_group), paste0(dir.results, "Simple_Effects_Sex_by_Group.csv"))
write.csv(as.data.frame(simple_effects_sex), paste0(dir.results, "Simple_Effects_Group_by_Sex.csv"))

# ============================================================================
# 3. Calculate Cohen's d effect sizes
# ============================================================================
cat("\n--- Cohen's d Effect Sizes ---\n")

library(effsize)

# Overall sex effect
cohens_d_sex <- cohen.d(model_data$pseudotime, model_data$sex)
cat("\nCohen's d for Sex (overall):\n")
print(cohens_d_sex)

# Overall group effect
cohens_d_group <- cohen.d(model_data$pseudotime, model_data$group_label)
cat("\nCohen's d for Diabetes Status (overall):\n")
print(cohens_d_group)

# Sex effect within Lean Control
lean_data <- model_data %>% filter(group_label == "Lean Control")
cohens_d_sex_lean <- cohen.d(lean_data$pseudotime, lean_data$sex)
cat("\nCohen's d for Sex in Lean Control:\n")
print(cohens_d_sex_lean)

# Sex effect within T2D
t2d_data <- model_data %>% filter(group_label == "T2D")
cohens_d_sex_t2d <- cohen.d(t2d_data$pseudotime, t2d_data$sex)
cat("\nCohen's d for Sex in T2D:\n")
print(cohens_d_sex_t2d)

# Group effect within Female
female_data <- model_data %>% filter(sex == "Female")
cohens_d_group_female <- cohen.d(female_data$pseudotime, female_data$group_label)
cat("\nCohen's d for Diabetes Status in Females:\n")
print(cohens_d_group_female)

# Group effect within Male
male_data <- model_data %>% filter(sex == "Male")
cohens_d_group_male <- cohen.d(male_data$pseudotime, male_data$group_label)
cat("\nCohen's d for Diabetes Status in Males:\n")
print(cohens_d_group_male)

# Compile effect sizes
effect_sizes_df <- data.frame(
  Comparison = c("Sex (Overall)", "Group (Overall)", 
                 "Sex in Lean Control", "Sex in T2D",
                 "Group in Females", "Group in Males"),
  Cohens_d = c(cohens_d_sex$estimate, cohens_d_group$estimate,
               cohens_d_sex_lean$estimate, cohens_d_sex_t2d$estimate,
               cohens_d_group_female$estimate, cohens_d_group_male$estimate),
  CI_lower = c(cohens_d_sex$conf.int[1], cohens_d_group$conf.int[1],
               cohens_d_sex_lean$conf.int[1], cohens_d_sex_t2d$conf.int[1],
               cohens_d_group_female$conf.int[1], cohens_d_group_male$conf.int[1]),
  CI_upper = c(cohens_d_sex$conf.int[2], cohens_d_group$conf.int[2],
               cohens_d_sex_lean$conf.int[2], cohens_d_sex_t2d$conf.int[2],
               cohens_d_group_female$conf.int[2], cohens_d_group_male$conf.int[2]),
  Magnitude = c(cohens_d_sex$magnitude, cohens_d_group$magnitude,
                cohens_d_sex_lean$magnitude, cohens_d_sex_t2d$magnitude,
                cohens_d_group_female$magnitude, cohens_d_group_male$magnitude)
)

print(effect_sizes_df)
write.csv(effect_sizes_df, paste0(dir.results, "Effect_Sizes_Cohens_d.csv"), row.names = FALSE)

# ============================================================================
# 4. Visualization: Interaction Plot with Error Bars
# ============================================================================
cat("\n--- Creating Interaction Plots ---\n")

# Calculate summary statistics for plotting
interaction_summary <- model_data %>%
  group_by(sex, group_label) %>%
  summarise(
    mean = mean(pseudotime),
    se = sd(pseudotime) / sqrt(n()),
    ci_lower = mean - 1.96 * se,
    ci_upper = mean + 1.96 * se,
    .groups = 'drop'
  )

# Interaction plot with lines
p_interaction1 <- ggplot(interaction_summary, aes(x = group_label, y = mean, 
                                                  color = sex, group = sex)) +
  geom_line(size = 1.5) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1, size = 1) +
  scale_color_manual(values = c("Female" = "#E64B35", "Male" = "#4DBBD5")) +
  theme_classic() +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  labs(title = "Interaction Plot: Sex × Diabetes Status on Pseudotime",
       subtitle = "Error bars represent 95% CI",
       x = "Diabetes Status",
       y = "Mean Pseudotime",
       color = "Sex")

print(p_interaction1)
ggsave(paste0(dir.results, "Interaction_Plot_Lines.pdf"), 
       plot = p_interaction1, width = 8, height = 6)

# Alternative: Bar plot with interaction
p_interaction2 <- ggplot(interaction_summary, aes(x = sex, y = mean, fill = group_label)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha = 0.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                position = position_dodge(0.9), width = 0.25, size = 0.8) +
  scale_fill_manual(values = c("T2D" = "#DC0000", "Lean Control" = "#00A087")) +
  theme_classic() +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  labs(title = "Mean Pseudotime by Sex and Diabetes Status",
       subtitle = "Error bars represent 95% CI",
       x = "Sex",
       y = "Mean Pseudotime",
       fill = "Diabetes Status")

print(p_interaction2)
ggsave(paste0(dir.results, "Interaction_Plot_Bars.pdf"), 
       plot = p_interaction2, width = 8, height = 6)

# ============================================================================
# 5. Violin Plot with Statistical Annotations (Enhanced)
# ============================================================================
cat("\n--- Creating Enhanced Violin Plot ---\n")

# Get emmeans for plotting
emm_df <- as.data.frame(emm)

p_violin_enhanced <- ggplot(model_data, aes(x = sex, y = pseudotime, fill = group_label)) +
  geom_violin(alpha = 0.6, trim = FALSE, position = position_dodge(0.8)) +
  geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA, 
               position = position_dodge(0.8)) +
  # Add mean points from emmeans
  geom_point(data = emm_df, aes(x = sex, y = emmean, group = group_label),
             position = position_dodge(0.8), size = 3, shape = 23, 
             fill = "white", color = "black") +
  scale_fill_manual(values = c("T2D" = "#DC0000", "Lean Control" = "#00A087")) +
  theme_classic() +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "right") +
  labs(title = "Pseudotime Distribution: Sex × Diabetes Status Interaction",
       subtitle = "Diamond points = estimated marginal means",
       x = "Sex",
       y = "Pseudotime",
       fill = "Diabetes Status")

print(p_violin_enhanced)
ggsave(paste0(dir.results, "Violin_Plot_Enhanced_with_Stats.pdf"), 
       plot = p_violin_enhanced, width = 10, height = 7)

# Alternative: Create separate violin plots with p-values for each comparison
# Plot A: Sex comparison within each diabetes group
p_violin_sex_by_group <- ggplot(model_data, aes(x = sex, y = pseudotime, fill = sex)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = c("Female" = "#E64B35", "Male" = "#4DBBD5")) +
  facet_wrap(~group_label) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_classic() +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none") +
  labs(title = "Sex Differences in Pseudotime within Each Diabetes Group",
       x = "Sex",
       y = "Pseudotime")

print(p_violin_sex_by_group)
ggsave(paste0(dir.results, "Violin_Plot_Sex_by_Group.pdf"), 
       plot = p_violin_sex_by_group, width = 10, height = 5)

# Plot B: Group comparison within each sex
p_violin_group_by_sex <- ggplot(model_data, aes(x = group_label, y = pseudotime, fill = group_label)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = c("T2D" = "#DC0000", "Lean Control" = "#00A087")) +
  facet_wrap(~sex) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_classic() +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none") +
  labs(title = "Diabetes Status Differences in Pseudotime within Each Sex",
       x = "Diabetes Status",
       y = "Pseudotime")

print(p_violin_group_by_sex)
ggsave(paste0(dir.results, "Violin_Plot_Group_by_Sex.pdf"), 
       plot = p_violin_group_by_sex, width = 10, height = 5)

# Plot C: All four groups side-by-side with manual annotation
# Calculate positions for significance bars
max_y <- max(model_data$pseudotime)
sig_height <- max_y * 1.1

p_violin_all_groups <- ggplot(model_data, aes(x = interaction(sex, group_label), 
                                              y = pseudotime, 
                                              fill = group_label)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = NA) +
  # Add mean points from emmeans
  geom_point(data = emm_df, 
             aes(x = interaction(sex, group_label), y = emmean),
             size = 3, shape = 23, fill = "white", color = "black") +
  scale_fill_manual(values = c("T2D" = "#DC0000", "Lean Control" = "#00A087")) +
  scale_x_discrete(labels = c("Female\nLean Control", "Male\nLean Control",
                              "Female\nT2D", "Male\nT2D")) +
  theme_classic() +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "right") +
  labs(title = "Pseudotime Distribution Across All Groups",
       subtitle = "Diamond points = estimated marginal means",
       x = "",
       y = "Pseudotime",
       fill = "Diabetes Status")

print(p_violin_all_groups)
ggsave(paste0(dir.results, "Violin_Plot_All_Groups.pdf"), 
       plot = p_violin_all_groups, width = 10, height = 7)

# ============================================================================
# 5b. Create a comprehensive comparison plot with manual p-value annotations
# ============================================================================

# Convert emmGrid objects to data frames
simple_effects_group_df <- as.data.frame(simple_effects_group)
simple_effects_sex_df <- as.data.frame(simple_effects_sex)

# Extract simple effects p-values
sex_in_lean <- simple_effects_group_df %>% 
  filter(group_label == "Lean Control") %>% 
  pull(p.value)

sex_in_t2d <- simple_effects_group_df %>% 
  filter(group_label == "T2D") %>% 
  pull(p.value)

group_in_female <- simple_effects_sex_df %>% 
  filter(sex == "Female") %>% 
  pull(p.value)

group_in_male <- simple_effects_sex_df %>% 
  filter(sex == "Male") %>% 
  pull(p.value)

# Format p-values
format_p <- function(p) {
  if (is.null(p) || length(p) == 0) return("p = NA")
  if (p < 0.001) return("p < 0.001***")
  else if (p < 0.01) return(sprintf("p = %.3f**", p))
  else if (p < 0.05) return(sprintf("p = %.3f*", p))
  else return(sprintf("p = %.3f", p))
}

# Get interaction p-value
interaction_p <- anova_results["sex:group_label", "Pr(>F)"]

# Create annotation data frame
annotations <- data.frame(
  label = c(
    paste0("Sex effect in Lean Control: ", format_p(sex_in_lean)),
    paste0("Sex effect in T2D: ", format_p(sex_in_t2d)),
    paste0("Group effect in Females: ", format_p(group_in_female)),
    paste0("Group effect in Males: ", format_p(group_in_male))
  ),
  x = 1,
  y = max_y * c(1.25, 1.2, 1.15, 1.1)
)

# Print the statistics for review
cat("\n=== Simple Effects Summary ===\n")
cat("\nSex effect in Lean Control: ", format_p(sex_in_lean), "\n")
cat("Sex effect in T2D: ", format_p(sex_in_t2d), "\n")
cat("Group effect in Females: ", format_p(group_in_female), "\n")
cat("Group effect in Males: ", format_p(group_in_male), "\n")
cat("Interaction (Sex × Group): ", format_p(interaction_p), "\n")

p_violin_annotated <- p_violin_all_groups +
  annotate("text", x = 1, y = max_y * 1.35, 
           label = paste0("Interaction (Sex × Diabetes): ", format_p(interaction_p)),
           hjust = 0, size = 5, fontface = "bold") +
  annotate("text", x = 1, y = annotations$y, 
           label = annotations$label,
           hjust = 0, size = 4)

print(p_violin_annotated)
ggsave(paste0(dir.results, "Violin_Plot_Fully_Annotated.pdf"), 
       plot = p_violin_annotated, width = 12, height = 8)

# ============================================================================
# 6. Model Diagnostics
# ============================================================================
cat("\n--- Model Diagnostics ---\n")

pdf(paste0(dir.results, "Model_Diagnostics.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2))
plot(lm_model)
dev.off()

# Check assumptions
n_residuals <- length(residuals(lm_model))
cat("\nNumber of observations:", n_residuals, "\n")

# Check normality of residuals
if (n_residuals <= 5000) {
  # Use Shapiro-Wilk test for smaller samples
  shapiro_test <- shapiro.test(residuals(lm_model))
  cat("\nShapiro-Wilk test for normality of residuals:\n")
  print(shapiro_test)
  
  # Save results
  sink(paste0(dir.results, "Normality_Test_Results.txt"))
  cat("Shapiro-Wilk Test for Normality of Residuals\n")
  cat("=============================================\n")
  print(shapiro_test)
  sink()
} else {
  # Use Kolmogorov-Smirnov test for larger samples
  cat("\nSample size exceeds 5000. Using Kolmogorov-Smirnov test instead.\n")
  ks_test <- ks.test(residuals(lm_model), "pnorm", 
                     mean = mean(residuals(lm_model)), 
                     sd = sd(residuals(lm_model)))
  cat("\nKolmogorov-Smirnov test for normality of residuals:\n")
  print(ks_test)
  
  # Also use Anderson-Darling test (from nortest package)
  if (requireNamespace("nortest", quietly = TRUE)) {
    library(nortest)
    ad_test <- ad.test(residuals(lm_model))
    cat("\nAnderson-Darling test for normality of residuals:\n")
    print(ad_test)
    
    # Save results
    sink(paste0(dir.results, "Normality_Test_Results.txt"))
    cat("Normality Tests for Residuals (n > 5000)\n")
    cat("==========================================\n\n")
    cat("Kolmogorov-Smirnov Test:\n")
    print(ks_test)
    cat("\n\nAnderson-Darling Test:\n")
    print(ad_test)
    sink()
  } else {
    # Save results
    sink(paste0(dir.results, "Normality_Test_Results.txt"))
    cat("Normality Test for Residuals (n > 5000)\n")
    cat("==========================================\n\n")
    cat("Kolmogorov-Smirnov Test:\n")
    print(ks_test)
    cat("\n\nNote: Install 'nortest' package for Anderson-Darling test: install.packages('nortest')\n")
    sink()
  }
  
  # Visual assessment is more important for large samples
  cat("\nNote: For large sample sizes (n > 5000), visual inspection of Q-Q plots")
  cat("\nis more informative than formal normality tests.\n")
}

# Create additional diagnostic plots
pdf(paste0(dir.results, "Extended_Model_Diagnostics.pdf"), width = 12, height = 10)

# Layout for 6 plots
par(mfrow = c(2, 3))

# 1. Residuals vs Fitted
plot(fitted(lm_model), residuals(lm_model),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

# 2. Q-Q plot
qqnorm(residuals(lm_model))
qqline(residuals(lm_model), col = "red")

# 3. Scale-Location
plot(fitted(lm_model), sqrt(abs(residuals(lm_model))),
     xlab = "Fitted Values", ylab = "√|Residuals|",
     main = "Scale-Location")
abline(h = mean(sqrt(abs(residuals(lm_model)))), col = "red", lty = 2)

# 4. Residuals vs Leverage
plot(lm_model, which = 5)

# 5. Histogram of residuals
hist(residuals(lm_model), breaks = 50, 
     main = "Histogram of Residuals",
     xlab = "Residuals", col = "lightblue")

# 6. Residuals by group
boxplot(residuals(lm_model) ~ model_data$sex:model_data$group_label,
        main = "Residuals by Group",
        xlab = "Group", ylab = "Residuals",
        las = 2, col = c("#00A087", "#DC0000", "#00A087", "#DC0000"))
abline(h = 0, col = "red", lty = 2)

dev.off()

# Levene's test for homogeneity of variance
library(car)
levene_test <- leveneTest(pseudotime ~ sex * group_label, data = model_data)
cat("\nLevene's test for homogeneity of variance:\n")
print(levene_test)

# Save Levene's test results
sink(paste0(dir.results, "Levene_Test_Results.txt"))
cat("Levene's Test for Homogeneity of Variance\n")
cat("==========================================\n")
print(levene_test)
sink()

# Additional diagnostics: Check for influential points
cat("\n--- Checking for Influential Observations ---\n")

# Cook's distance
cooks_d <- cooks.distance(lm_model)
influential_threshold <- 4 / n_residuals
n_influential <- sum(cooks_d > influential_threshold)
cat("\nNumber of potentially influential observations (Cook's D > 4/n):", n_influential, "\n")

if (n_influential > 0 && n_influential < 100) {
  cat("Top 10 influential observations by Cook's distance:\n")
  top_influential <- sort(cooks_d, decreasing = TRUE)[1:min(10, n_influential)]
  print(top_influential)
}

# Plot Cook's distance
pdf(paste0(dir.results, "Cooks_Distance.pdf"), width = 10, height = 6)
plot(cooks_d, type = "h", main = "Cook's Distance",
     ylab = "Cook's Distance", xlab = "Observation Index")
abline(h = influential_threshold, col = "red", lty = 2)
legend("topright", legend = paste("Threshold: 4/n =", round(influential_threshold, 4)),
       col = "red", lty = 2)
dev.off()

# VIF for multicollinearity (if applicable)
if (requireNamespace("car", quietly = TRUE)) {
  # VIF can only be calculated if we have multiple predictors
  # For interaction models, we check VIF on main effects model
  lm_main <- lm(pseudotime ~ sex + group_label, data = model_data)
  vif_values <- car::vif(lm_main)
  cat("\nVariance Inflation Factors (VIF) for main effects:\n")
  print(vif_values)
  
  sink(paste0(dir.results, "VIF_Results.txt"))
  cat("Variance Inflation Factors\n")
  cat("==========================\n")
  print(vif_values)
  cat("\nNote: VIF > 10 indicates potential multicollinearity problems\n")
  sink()
}

# Summary of diagnostic checks
cat("\n=== Diagnostic Summary ===\n")
cat("1. Residual plots saved to: Model_Diagnostics.pdf and Extended_Model_Diagnostics.pdf\n")
cat("2. Normality test results saved to: Normality_Test_Results.txt\n")
cat("3. Homogeneity of variance test saved to: Levene_Test_Results.txt\n")
cat("4. Cook's distance plot saved to: Cooks_Distance.pdf\n")
cat("5. Number of potentially influential observations:", n_influential, "\n")


# ============================================================================
# 7. Summary Table for Publication
# ============================================================================
cat("\n--- Creating Summary Table ---\n")

summary_table <- model_data %>%
  group_by(sex, group_label) %>%
  summarise(
    N = n(),
    Mean = mean(pseudotime),
    SD = sd(pseudotime),
    Median = median(pseudotime),
    IQR = IQR(pseudotime),
    Min = min(pseudotime),
    Max = max(pseudotime),
    .groups = 'drop'
  ) %>%
  arrange(group_label, sex)

print(summary_table)
write.csv(summary_table, paste0(dir.results, "Summary_Table_Publication.csv"), 
          row.names = FALSE)

# Create formatted table
library(knitr)
library(kableExtra)

formatted_table <- summary_table %>%
  kable(digits = 3, caption = "Pseudotime Summary Statistics by Sex and Diabetes Status") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))

cat(formatted_table, file = paste0(dir.results, "Summary_Table_Formatted.html"))

# ============================================================================
# 8. Create Combined Results Plot
# ============================================================================
cat("\n--- Creating Combined Results Plot ---\n")

combined_stats_plot <- (p_interaction1 | p_interaction2) / 
  (p_violin_sex_by_group | p_violin_group_by_sex)

print(combined_stats_plot)
ggsave(paste0(dir.results, "Combined_Statistical_Analysis.pdf"), 
       plot = combined_stats_plot, width = 16, height = 12)

# Alternative combined plot with annotated version
combined_stats_plot2 <- (p_interaction1 | p_interaction2) / 
  p_violin_annotated

print(combined_stats_plot2)
ggsave(paste0(dir.results, "Combined_Statistical_Analysis_Annotated.pdf"), 
       plot = combined_stats_plot2, width = 16, height = 12)

cat("\n=== Statistical Analysis Complete ===\n")
cat("All results saved to:", dir.results, "\n")











# ============================================================================
# Pseudotime Density Plots Faceted by PT Cell Subtype
# ============================================================================
cat("\n=== Creating PT Subtype-Specific Density Plots ===\n")

# Check available PT subtypes
cat("\nAvailable PT cell subtypes:\n")
print(table(plot_df_clean$celltype))

# ============================================================================
# Plot 1: Density by Sex, faceted by PT subtype
# ============================================================================
cat("\nCreating density plot by sex, faceted by PT subtype...\n")

p_density_sex_by_subtype <- ggplot(plot_df_clean, aes(x = pseudotime, fill = sex)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("Female" = "#E64B35", "Male" = "#4DBBD5")) +
  facet_wrap(~celltype, scales = "free_y", ncol = 2) +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom") +
  labs(title = "Pseudotime Distribution by Sex Across PT Cell Subtypes",
       x = "Pseudotime",
       y = "Density",
       fill = "Sex")

print(p_density_sex_by_subtype)
ggsave(paste0(dir.results, "Density_Sex_by_PT_Subtype.pdf"), 
       plot = p_density_sex_by_subtype, width = 12, height = 10)

# ============================================================================
# Plot 2: Density by Diabetes Status, faceted by PT subtype
# ============================================================================
cat("\nCreating density plot by diabetes status, faceted by PT subtype...\n")

p_density_group_by_subtype <- ggplot(plot_df_clean, aes(x = pseudotime, fill = group_label)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("T2D" = "#DC0000", "Lean Control" = "#00A087")) +
  facet_wrap(~celltype, scales = "free_y", ncol = 2) +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom") +
  labs(title = "Pseudotime Distribution by Diabetes Status Across PT Cell Subtypes",
       x = "Pseudotime",
       y = "Density",
       fill = "Diabetes Status")

print(p_density_group_by_subtype)
ggsave(paste0(dir.results, "Density_Group_by_PT_Subtype.pdf"), 
       plot = p_density_group_by_subtype, width = 12, height = 10)

# ============================================================================
# Plot 3: Density by Sex within each Diabetes Status, faceted by PT subtype
# ============================================================================
cat("\nCreating density plot by sex within diabetes status, faceted by PT subtype...\n")

p_density_sex_group_by_subtype <- ggplot(plot_df_clean, aes(x = pseudotime, fill = sex)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("Female" = "#E64B35", "Male" = "#4DBBD5")) +
  facet_grid(celltype ~ group_label, scales = "free_y") +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom") +
  labs(title = "Pseudotime Distribution by Sex and Diabetes Status Across PT Cell Subtypes",
       x = "Pseudotime",
       y = "Density",
       fill = "Sex")

print(p_density_sex_group_by_subtype)
ggsave(paste0(dir.results, "Density_Sex_and_Group_by_PT_Subtype.pdf"), 
       plot = p_density_sex_group_by_subtype, width = 12, height = 14)

# ============================================================================
# Plot 4: Ridgeline plot - PT subtypes by Sex
# ============================================================================
cat("\nCreating ridgeline plot of PT subtypes by sex...\n")

p_ridge_subtype_by_sex <- ggplot(plot_df_clean, aes(x = pseudotime, y = celltype, fill = sex)) +
  geom_density_ridges(alpha = 0.7, scale = 2.5, rel_min_height = 0.01) +
  scale_fill_manual(values = c("Female" = "#E64B35", "Male" = "#4DBBD5")) +
  theme_ridges() +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 11),
        legend.position = "right") +
  labs(title = "PT Cell Subtype Distribution Along Pseudotime by Sex",
       x = "Pseudotime",
       y = "PT Cell Subtype",
       fill = "Sex")

print(p_ridge_subtype_by_sex)
ggsave(paste0(dir.results, "Ridgeline_PT_Subtype_by_Sex.pdf"), 
       plot = p_ridge_subtype_by_sex, width = 10, height = 8)

# ============================================================================
# Plot 5: Ridgeline plot - PT subtypes by Diabetes Status
# ============================================================================
cat("\nCreating ridgeline plot of PT subtypes by diabetes status...\n")

p_ridge_subtype_by_group <- ggplot(plot_df_clean, aes(x = pseudotime, y = celltype, fill = group_label)) +
  geom_density_ridges(alpha = 0.7, scale = 2.5, rel_min_height = 0.01) +
  scale_fill_manual(values = c("T2D" = "#DC0000", "Lean Control" = "#00A087")) +
  theme_ridges() +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 11),
        legend.position = "right") +
  labs(title = "PT Cell Subtype Distribution Along Pseudotime by Diabetes Status",
       x = "Pseudotime",
       y = "PT Cell Subtype",
       fill = "Diabetes Status")

print(p_ridge_subtype_by_group)
ggsave(paste0(dir.results, "Ridgeline_PT_Subtype_by_Group.pdf"), 
       plot = p_ridge_subtype_by_group, width = 10, height = 8)

# ============================================================================
# Plot 6: Faceted ridgeline - PT subtypes by Sex and Diabetes Status
# ============================================================================
cat("\nCreating faceted ridgeline plot of PT subtypes...\n")

p_ridge_subtype_faceted <- ggplot(plot_df_clean, aes(x = pseudotime, y = celltype, 
                                                     fill = interaction(sex, group_label))) +
  geom_density_ridges(alpha = 0.7, scale = 2, rel_min_height = 0.01) +
  scale_fill_manual(values = c("Female.Lean Control" = "#E64B35",
                               "Male.Lean Control" = "#4DBBD5",
                               "Female.T2D" = "#DC0000",
                               "Male.T2D" = "#00A087"),
                    labels = c("Female Lean Control", "Male Lean Control",
                               "Female T2D", "Male T2D")) +
  theme_ridges() +
  theme(text = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_blank()) +
  labs(title = "PT Cell Subtype Distribution by Sex and Diabetes Status",
       x = "Pseudotime",
       y = "PT Cell Subtype")

print(p_ridge_subtype_faceted)
ggsave(paste0(dir.results, "Ridgeline_PT_Subtype_All_Groups.pdf"), 
       plot = p_ridge_subtype_faceted, width = 11, height = 9)

# ============================================================================
# Plot 7: Heatmap of mean pseudotime by PT subtype and group
# ============================================================================
cat("\nCreating heatmap of mean pseudotime...\n")

# Calculate mean pseudotime for each combination
heatmap_data <- plot_df_clean %>%
  group_by(celltype, sex, group_label) %>%
  summarise(
    mean_pseudotime = mean(pseudotime),
    n_cells = n(),
    .groups = 'drop'
  ) %>%
  mutate(group_combo = paste(sex, group_label, sep = "\n"))

# Create heatmap
p_heatmap <- ggplot(heatmap_data, aes(x = group_combo, y = celltype, fill = mean_pseudotime)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = round(mean_pseudotime, 1)), color = "white", size = 4, fontface = "bold") +
  scale_fill_viridis_c(option = "plasma") +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.grid = element_blank(),
        legend.position = "right") +
  labs(title = "Mean Pseudotime Across PT Cell Subtypes by Sex and Diabetes Status",
       x = "",
       y = "PT Cell Subtype",
       fill = "Mean\nPseudotime")

print(p_heatmap)
ggsave(paste0(dir.results, "Heatmap_Mean_Pseudotime_PT_Subtype.pdf"), 
       plot = p_heatmap, width = 10, height = 8)

# Save the summary data
write.csv(heatmap_data, paste0(dir.results, "Mean_Pseudotime_by_PT_Subtype.csv"), 
          row.names = FALSE)

# ============================================================================
# Plot 8: Violin plots by PT subtype
# ============================================================================
cat("\nCreating violin plots by PT subtype...\n")

p_violin_by_subtype <- ggplot(plot_df_clean, aes(x = celltype, y = pseudotime, 
                                                 fill = interaction(sex, group_label))) +
  geom_violin(alpha = 0.7, position = position_dodge(0.8), trim = FALSE) +
  geom_boxplot(width = 0.15, position = position_dodge(0.8), 
               alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = c("Female.Lean Control" = "#E64B35",
                               "Male.Lean Control" = "#4DBBD5",
                               "Female.T2D" = "#DC0000",
                               "Male.T2D" = "#00A087"),
                    labels = c("Female Lean Control", "Male Lean Control",
                               "Female T2D", "Male T2D")) +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        legend.position = "bottom",
        legend.title = element_blank()) +
  labs(title = "Pseudotime Distribution Across PT Cell Subtypes",
       x = "PT Cell Subtype",
       y = "Pseudotime")

print(p_violin_by_subtype)
ggsave(paste0(dir.results, "Violin_Pseudotime_by_PT_Subtype.pdf"), 
       plot = p_violin_by_subtype, width = 12, height = 8)

# ============================================================================
# Plot 9: Cell count distribution across PT subtypes
# ============================================================================
cat("\nCreating cell count distribution plot...\n")

cell_counts <- plot_df_clean %>%
  group_by(celltype, sex, group_label) %>%
  summarise(n_cells = n(), .groups = 'drop')

p_cell_counts <- ggplot(cell_counts, aes(x = celltype, y = n_cells, 
                                         fill = interaction(sex, group_label))) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("Female.Lean Control" = "#E64B35",
                               "Male.Lean Control" = "#4DBBD5",
                               "Female.T2D" = "#DC0000",
                               "Male.T2D" = "#00A087"),
                    labels = c("Female Lean Control", "Male Lean Control",
                               "Female T2D", "Male T2D")) +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        legend.position = "bottom",
        legend.title = element_blank()) +
  labs(title = "Cell Count Distribution Across PT Cell Subtypes",
       x = "PT Cell Subtype",
       y = "Number of Cells")

print(p_cell_counts)
ggsave(paste0(dir.results, "Cell_Counts_by_PT_Subtype.pdf"), 
       plot = p_cell_counts, width = 12, height = 7)

# Save cell count data
write.csv(cell_counts, paste0(dir.results, "Cell_Counts_by_PT_Subtype.csv"), 
          row.names = FALSE)

# ============================================================================
# Plot 10: Combined panel figure for PT subtypes
# ============================================================================
cat("\nCreating combined panel figure...\n")

combined_subtype_plot <- (p_density_sex_by_subtype | p_density_group_by_subtype) /
  (p_ridge_subtype_by_sex | p_ridge_subtype_by_group)

print(combined_subtype_plot)
ggsave(paste0(dir.results, "Combined_PT_Subtype_Analysis.pdf"), 
       plot = combined_subtype_plot, width = 18, height = 14)

# ============================================================================
# Statistical analysis by PT subtype
# ============================================================================
cat("\n--- Statistical Analysis by PT Subtype ---\n")

# Test for sex and group differences within each PT subtype
subtype_stats <- plot_df_clean %>%
  group_by(celltype) %>%
  summarise(
    n_cells = n(),
    # Sex effect (Wilcoxon test)
    sex_p_value = wilcox.test(pseudotime ~ sex)$p.value,
    # Group effect (Wilcoxon test)
    group_p_value = wilcox.test(pseudotime ~ group_label)$p.value,
    # Mean pseudotime by sex
    mean_female = mean(pseudotime[sex == "Female"]),
    mean_male = mean(pseudotime[sex == "Male"]),
    # Mean pseudotime by group
    mean_lean = mean(pseudotime[group_label == "Lean Control"]),
    mean_t2d = mean(pseudotime[group_label == "T2D"]),
    .groups = 'drop'
  ) %>%
  mutate(
    sex_p_adj = p.adjust(sex_p_value, method = "bonferroni"),
    group_p_adj = p.adjust(group_p_value, method = "bonferroni"),
    sex_significant = sex_p_adj < 0.05,
    group_significant = group_p_adj < 0.05
  )

cat("\nStatistical tests by PT subtype:\n")
print(subtype_stats)

write.csv(subtype_stats, paste0(dir.results, "Statistical_Tests_by_PT_Subtype.csv"), 
          row.names = FALSE)

# ============================================================================
# Summary statistics table by PT subtype
# ============================================================================
summary_by_subtype <- plot_df_clean %>%
  group_by(celltype, sex, group_label) %>%
  summarise(
    n_cells = n(),
    mean_pseudotime = mean(pseudotime),
    sd_pseudotime = sd(pseudotime),
    median_pseudotime = median(pseudotime),
    q25 = quantile(pseudotime, 0.25),
    q75 = quantile(pseudotime, 0.75),
    .groups = 'drop'
  ) %>%
  arrange(celltype, group_label, sex)

cat("\nSummary statistics by PT subtype:\n")
print(summary_by_subtype)

write.csv(summary_by_subtype, paste0(dir.results, "Summary_Statistics_by_PT_Subtype.csv"), 
          row.names = FALSE)

cat("\n=== PT Subtype Analysis Complete ===\n")
cat("All plots and statistics saved to:", dir.results, "\n")









