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
































####Overall comparisons

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Define directories
base_dir <- "C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis"
lc_dir <- file.path(base_dir, "LeanControl_Only")
t2d_dir <- file.path(base_dir, "T2D_Only")
results_dir <- file.path(base_dir, "Sex_Comparison_Results")

# Create results directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
  cat(sprintf("Created results directory: %s\n\n", results_dir))
} else {
  cat(sprintf("Using existing results directory: %s\n\n", results_dir))
}

# Get all cell type files from each directory
lc_files <- list.files(lc_dir, pattern = "Full_NEBULA_.*_LC_pooledoffset\\.csv", full.names = TRUE)
t2d_files <- list.files(t2d_dir, pattern = "Full_NEBULA_.*_T2D_pooledoffset\\.csv", full.names = TRUE)

cat(sprintf("Found %d LC files\n", length(lc_files)))
cat(sprintf("Found %d T2D files\n\n", length(t2d_files)))

# Extract cell types from filenames
extract_celltype <- function(filepath) {
  filename <- basename(filepath)
  # Extract cell type between "Full_NEBULA_" and "_LC" or "_T2D"
  gsub("Full_NEBULA_(.+?)_(LC|T2D)_pooledoffset\\.csv", "\\1", filename)
}

# Load all data
load_and_annotate <- function(files, condition) {
  map_dfr(files, function(f) {
    cat(sprintf("Loading: %s\n", basename(f)))
    read_csv(f, show_col_types = FALSE) %>%
      mutate(
        celltype = extract_celltype(f),
        condition = condition
      )
  })
}

cat("Loading Lean Control files...\n")
lc_data <- load_and_annotate(lc_files, "LC")

cat("\nLoading T2D files...\n")
t2d_data <- load_and_annotate(t2d_files, "T2D")

cat(sprintf("\nTotal LC observations: %d\n", nrow(lc_data)))
cat(sprintf("Total T2D observations: %d\n", nrow(t2d_data)))
cat(sprintf("Cell types in LC: %s\n", paste(unique(lc_data$celltype), collapse = ", ")))
cat(sprintf("Cell types in T2D: %s\n\n", paste(unique(t2d_data$celltype), collapse = ", ")))

# ============================================================================
# 2. PREPARE SEX EFFECT DATA
# ============================================================================

# Focus on sex effect (sexMale coefficient)
# Positive logFC = higher in males, Negative logFC = higher in females (protective)

prepare_sex_effects <- function(df) {
  df %>%
    dplyr::select(
      summary.gene,
      celltype,
      condition,
      logFC_sex = summary.logFC_sexMale,
      se_sex = summary.se_sexMale,
      p_sex = summary.p_sexMale
    ) %>%
    mutate(
      padj_sex = p.adjust(p_sex, method = "BH"),
      sig = padj_sex < 0.05,
      direction = case_when(
        logFC_sex > 0 ~ "Male-biased",
        logFC_sex < 0 ~ "Female-biased",
        TRUE ~ "No difference"
      )
    )
}

lc_sex <- prepare_sex_effects(lc_data)
t2d_sex <- prepare_sex_effects(t2d_data)

# Combine for comparison
combined <- bind_rows(lc_sex, t2d_sex) %>%
  pivot_wider(
    id_cols = c(summary.gene, celltype),
    names_from = condition,
    values_from = c(logFC_sex, se_sex, p_sex, padj_sex, sig, direction)
  )

# ============================================================================
# 3. CATEGORIZE GENES BY SEX EFFECT PATTERNS
# ============================================================================

combined <- combined %>%
  mutate(
    pattern = case_when(
      # Female-protective in both conditions
      sig_LC & sig_T2D & logFC_sex_LC < 0 & logFC_sex_T2D < 0 ~ 
        "Female-protective (both)",
      
      # Female-protective only in LC (lost in T2D)
      sig_LC & !sig_T2D & logFC_sex_LC < 0 ~ 
        "Female-protective lost in T2D",
      
      # Female-protective only in T2D (gained)
      !sig_LC & sig_T2D & logFC_sex_T2D < 0 ~ 
        "Female-protective gained in T2D",
      
      # Attenuated protection (both sig, both negative, but smaller in T2D)
      sig_LC & sig_T2D & logFC_sex_LC < 0 & logFC_sex_T2D < 0 & 
        abs(logFC_sex_T2D) < abs(logFC_sex_LC) ~ 
        "Attenuated female protection",
      
      # Male-biased in both
      sig_LC & sig_T2D & logFC_sex_LC > 0 & logFC_sex_T2D > 0 ~ 
        "Male-biased (both)",
      
      # Other patterns
      TRUE ~ "Other/Not significant"
    ),
    
    # Calculate delta (change in sex effect from LC to T2D)
    delta_logFC = logFC_sex_T2D - logFC_sex_LC,
    
    # Is the female protection attenuated? (becoming less negative)
    protection_attenuated = logFC_sex_LC < 0 & logFC_sex_T2D > logFC_sex_LC
  )

# ============================================================================
# 4. SUMMARY STATISTICS BY CELL TYPE
# ============================================================================

summary_by_celltype <- combined %>%
  group_by(celltype) %>%
  summarise(
    n_genes = n(),
    
    # Female-protective genes
    n_female_prot_LC = sum(sig_LC & logFC_sex_LC < 0, na.rm = TRUE),
    n_female_prot_T2D = sum(sig_T2D & logFC_sex_T2D < 0, na.rm = TRUE),
    pct_lost = 100 * sum(pattern == "Female-protective lost in T2D", na.rm = TRUE) / n_genes,
    
    # Attenuation metrics
    n_attenuated = sum(protection_attenuated, na.rm = TRUE),
    mean_attenuation = mean(delta_logFC[protection_attenuated], na.rm = TRUE),
    
    # Male-biased genes
    n_male_LC = sum(sig_LC & logFC_sex_LC > 0, na.rm = TRUE),
    n_male_T2D = sum(sig_T2D & logFC_sex_T2D > 0, na.rm = TRUE)
  ) %>%
  arrange(desc(n_female_prot_LC))

print("Summary by cell type:")
print(summary_by_celltype)

# ============================================================================
# 5. VISUALIZATIONS
# ============================================================================

# 5a. Scatter plot: LC vs T2D sex effects by cell type
for (ct in unique(combined$celltype)) {
  p <- combined %>%
    filter(celltype == ct) %>%
    filter(logFC_sex_LC >= -15 & logFC_sex_LC <= 15) %>%
    filter(logFC_sex_T2D >= -15 & logFC_sex_T2D <= 15) %>%
    ggplot(aes(x = logFC_sex_LC, y = logFC_sex_T2D)) +
    geom_point(aes(color = pattern, alpha = sig_LC | sig_T2D), size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    scale_alpha_manual(values = c(0.3, 0.8)) +
    labs(
      title = paste("Sex Effects:", ct),
      subtitle = "Comparing Lean Controls vs T2D",
      x = "Log2FC (Male vs Female) in Lean Controls",
      y = "Log2FC (Male vs Female) in T2D",
      color = "Pattern",
      caption = "Points above diagonal = increased male bias in T2D\nPoints below diagonal = increased female protection in T2D"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p)
  
  # Save plot
  ggsave(file.path(results_dir, paste0("sex_comparison_", ct, ".png")), 
         p, width = 10, height = 8)
}

# 5b. Heatmap of delta (change in sex effect)
top_attenuated <- combined %>%
  filter(protection_attenuated) %>%
  group_by(celltype) %>%
  slice_max(abs(delta_logFC), n = 20) %>%
  ungroup()

if (nrow(top_attenuated) > 0) {
  heatmap_data <- top_attenuated %>%
    dplyr::select(summary.gene, celltype, delta_logFC) %>%
    pivot_wider(names_from = celltype, values_from = delta_logFC) %>%
    column_to_rownames("summary.gene")
  
  # Remove rows or columns with all NAs
  heatmap_data <- heatmap_data[rowSums(is.na(heatmap_data)) < ncol(heatmap_data), ]
  heatmap_data <- heatmap_data[, colSums(is.na(heatmap_data)) < nrow(heatmap_data)]
  
  heatmap_data <- as.matrix(heatmap_data)
  # Check if we still have data
  if (nrow(heatmap_data) > 1 && ncol(heatmap_data) > 1) {
    # Replace any remaining NAs with 0 for visualization
    heatmap_data[is.na(heatmap_data)] <- 0
    
    # Remove any infinite values
    heatmap_data[is.infinite(heatmap_data)] <- 0
    
    pheatmap(
      heatmap_data,
      main = "Top Genes with Attenuated Female Protection in T2D",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      color = colorRampPalette(c("blue", "white", "red"))(50),
      breaks = seq(-max(abs(heatmap_data), na.rm = TRUE), 
                   max(abs(heatmap_data), na.rm = TRUE), 
                   length.out = 51),
      na_col = "grey"
    )
  } else {
    cat("Not enough data for heatmap after filtering NAs\n")
  }
}

# 5c. Bar plot: Pattern distribution by cell type
pattern_summary <- combined %>%
  count(celltype, pattern) %>%
  filter(pattern != "Other/Not significant") %>%
  ggplot(aes(x = celltype, y = n, fill = pattern)) +
  geom_col(position = "dodge") +
  labs(
    title = "Sex Effect Patterns Across Cell Types",
    x = "Cell Type",
    y = "Number of Genes",
    fill = "Pattern"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(pattern_summary)
ggsave(file.path(results_dir, "pattern_summary_by_celltype.png"), 
       pattern_summary, width = 12, height = 6)

# ============================================================================
# 6. GO PATHWAY ANALYSIS - CELL TYPE BY CELL TYPE
# ============================================================================

cat("\n\n=== RUNNING GO PATHWAY ANALYSIS ===\n\n")

# Create GO results directory
go_dir <- file.path(results_dir, "GO_Pathway_Analysis")
if (!dir.exists(go_dir)) {
  dir.create(go_dir)
}

# Function to run GO enrichment
run_go_enrichment <- function(genes, analysis_name, universe = NULL) {
  if (length(genes) < 5) {
    cat(sprintf("  Skipping %s - too few genes (%d)\n", analysis_name, length(genes)))
    return(NULL)
  }
  
  tryCatch({
    ego <- enrichGO(
      gene = genes,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",  # Biological Process
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      universe = universe,
      readable = TRUE
    )
    
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      cat(sprintf("  ✓ %s: %d significant pathways\n", analysis_name, nrow(as.data.frame(ego))))
      return(ego)
    } else {
      cat(sprintf("  - %s: No significant pathways\n", analysis_name))
      return(NULL)
    }
  }, error = function(e) {
    cat(sprintf("  ✗ %s: Error - %s\n", analysis_name, e$message))
    return(NULL)
  })
}

# Analyze each cell type
all_go_results <- list()
go_comparison_data <- list()

for (ct in unique(combined$celltype)) {
  cat(sprintf("\n--- Analyzing %s ---\n", ct))
  
  ct_data <- combined %>% filter(celltype == ct)
  
  # Get universe (all genes tested in this cell type)
  universe_genes <- ct_data$summary.gene
  
  # Define gene sets for this cell type
  # 1. Female-protective genes in LC
  lc_female_prot <- ct_data %>%
    filter(sig_LC & logFC_sex_LC < 0) %>%
    pull(summary.gene)
  
  # 2. Female-protective genes in T2D
  t2d_female_prot <- ct_data %>%
    filter(sig_T2D & logFC_sex_T2D < 0) %>%
    pull(summary.gene)
  
  # 3. Female-protective genes LOST in T2D
  lost_protection <- ct_data %>%
    filter(pattern == "Female-protective lost in T2D") %>%
    pull(summary.gene)
  
  # 4. Female-protective genes with ATTENUATED protection
  attenuated_protection <- ct_data %>%
    filter(protection_attenuated & sig_LC) %>%
    pull(summary.gene)
  
  # 5. Male-biased genes in LC
  lc_male_bias <- ct_data %>%
    filter(sig_LC & logFC_sex_LC > 0) %>%
    pull(summary.gene)
  
  # 6. Male-biased genes in T2D
  t2d_male_bias <- ct_data %>%
    filter(sig_T2D & logFC_sex_T2D > 0) %>%
    pull(summary.gene)
  
  # Run GO enrichment for each gene set
  go_lc_female <- run_go_enrichment(lc_female_prot, 
                                    paste0(ct, " - Female-protective (LC)"), 
                                    universe_genes)
  
  go_t2d_female <- run_go_enrichment(t2d_female_prot, 
                                     paste0(ct, " - Female-protective (T2D)"), 
                                     universe_genes)
  
  go_lost <- run_go_enrichment(lost_protection, 
                               paste0(ct, " - Lost protection in T2D"), 
                               universe_genes)
  
  go_attenuated <- run_go_enrichment(attenuated_protection, 
                                     paste0(ct, " - Attenuated protection"), 
                                     universe_genes)
  
  go_lc_male <- run_go_enrichment(lc_male_bias, 
                                  paste0(ct, " - Male-biased (LC)"), 
                                  universe_genes)
  
  go_t2d_male <- run_go_enrichment(t2d_male_bias, 
                                   paste0(ct, " - Male-biased (T2D)"), 
                                   universe_genes)
  
  # Store results
  all_go_results[[ct]] <- list(
    lc_female = go_lc_female,
    t2d_female = go_t2d_female,
    lost = go_lost,
    attenuated = go_attenuated,
    lc_male = go_lc_male,
    t2d_male = go_t2d_male
  )
  
  # Create visualizations for this cell type if we have results
  ct_go_dir <- file.path(go_dir, ct)
  if (!dir.exists(ct_go_dir)) {
    dir.create(ct_go_dir)
  }
  
  # Plot LC female-protective pathways
  if (!is.null(go_lc_female) && nrow(as.data.frame(go_lc_female)) > 0) {
    p <- dotplot(go_lc_female, showCategory = 20, title = paste(ct, "- Female-protective Pathways (LC)"))
    ggsave(file.path(ct_go_dir, "GO_LC_female_protective.png"), p, width = 12, height = 8)
  }
  
  # Plot T2D female-protective pathways
  if (!is.null(go_t2d_female) && nrow(as.data.frame(go_t2d_female)) > 0) {
    p <- dotplot(go_t2d_female, showCategory = 20, title = paste(ct, "- Female-protective Pathways (T2D)"))
    ggsave(file.path(ct_go_dir, "GO_T2D_female_protective.png"), p, width = 12, height = 8)
  }
  
  # Plot lost protection pathways
  if (!is.null(go_lost) && nrow(as.data.frame(go_lost)) > 0) {
    p <- dotplot(go_lost, showCategory = 20, title = paste(ct, "- Pathways Lost in T2D"))
    ggsave(file.path(ct_go_dir, "GO_lost_protection.png"), p, width = 12, height = 8)
  }
  
  # Plot attenuated protection pathways
  if (!is.null(go_attenuated) && nrow(as.data.frame(go_attenuated)) > 0) {
    p <- dotplot(go_attenuated, showCategory = 20, title = paste(ct, "- Attenuated Protection Pathways"))
    ggsave(file.path(ct_go_dir, "GO_attenuated_protection.png"), p, width = 12, height = 8)
  }
  
  # Compare LC vs T2D female-protective pathways
  if (!is.null(go_lc_female) && !is.null(go_t2d_female)) {
    lc_pathways <- as.data.frame(go_lc_female) %>%
      dplyr::select(ID, Description, p.adjust, Count) %>%
      mutate(condition = "LC")
    
    t2d_pathways <- as.data.frame(go_t2d_female) %>%
      dplyr::select(ID, Description, p.adjust, Count) %>%
      mutate(condition = "T2D")
    
    # Combine and identify shared vs unique pathways
    comparison <- bind_rows(lc_pathways, t2d_pathways) %>%
      group_by(ID, Description) %>%
      summarise(
        in_LC = "LC" %in% condition,
        in_T2D = "T2D" %in% condition,
        lc_padj = ifelse(any(condition == "LC"), p.adjust[condition == "LC"][1], NA),
        t2d_padj = ifelse(any(condition == "T2D"), p.adjust[condition == "T2D"][1], NA),
        lc_count = ifelse(any(condition == "LC"), Count[condition == "LC"][1], NA),
        t2d_count = ifelse(any(condition == "T2D"), Count[condition == "T2D"][1], NA),
        .groups = "drop"
      ) %>%
      mutate(
        pathway_pattern = case_when(
          in_LC & in_T2D ~ "Shared",
          in_LC & !in_T2D ~ "LC only",
          !in_LC & in_T2D ~ "T2D only"
        )
      )
    
    go_comparison_data[[ct]] <- comparison
    
    # Save comparison
    write_csv(comparison, file.path(ct_go_dir, "pathway_comparison_LC_vs_T2D.csv"))
    
    # Visualize comparison
    if (nrow(comparison) > 0) {
      p_comp <- comparison %>%
        mutate(Description = str_wrap(Description, 50)) %>%
        head(30) %>%
        ggplot(aes(x = reorder(Description, lc_padj), 
                   y = -log10(lc_padj),
                   fill = pathway_pattern)) +
        geom_col() +
        coord_flip() +
        labs(
          title = paste(ct, "- Female-protective Pathway Comparison"),
          subtitle = "LC vs T2D",
          x = "Pathway",
          y = "-log10(adjusted p-value) in LC",
          fill = "Pattern"
        ) +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 8))
      
      ggsave(file.path(ct_go_dir, "pathway_comparison_barplot.png"), 
             p_comp, width = 12, height = 10)
    }
  }
  
  # Save all GO results as CSV for this cell type
  for (analysis_name in names(all_go_results[[ct]])) {
    go_result <- all_go_results[[ct]][[analysis_name]]
    if (!is.null(go_result) && nrow(as.data.frame(go_result)) > 0) {
      write_csv(
        as.data.frame(go_result),
        file.path(ct_go_dir, paste0("GO_", analysis_name, ".csv"))
      )
    }
  }
}

# ============================================================================
# 7. CROSS-CELL TYPE GO COMPARISON
# ============================================================================

cat("\n\n=== CROSS-CELL TYPE PATHWAY COMPARISON ===\n\n")

# Collect top pathways lost/attenuated across all cell types
all_lost_pathways <- map_dfr(names(all_go_results), function(ct) {
  if (!is.null(all_go_results[[ct]]$lost)) {
    as.data.frame(all_go_results[[ct]]$lost) %>%
      mutate(celltype = ct, category = "Lost in T2D") %>%
      head(10)
  }
})

all_attenuated_pathways <- map_dfr(names(all_go_results), function(ct) {
  if (!is.null(all_go_results[[ct]]$attenuated)) {
    as.data.frame(all_go_results[[ct]]$attenuated) %>%
      mutate(celltype = ct, category = "Attenuated") %>%
      head(10)
  }
})

# Combine and save
all_key_pathways <- bind_rows(all_lost_pathways, all_attenuated_pathways)

if (nrow(all_key_pathways) > 0) {
  write_csv(all_key_pathways, file.path(go_dir, "key_pathways_across_celltypes.csv"))
  
  # Visualize pathways appearing in multiple cell types
  pathway_counts <- all_key_pathways %>%
    count(Description, category, sort = TRUE) %>%
    filter(n > 1) %>%
    head(20)
  
  if (nrow(pathway_counts) > 0) {
    p_shared <- ggplot(pathway_counts, 
                       aes(x = reorder(str_wrap(Description, 40), n), 
                           y = n, 
                           fill = category)) +
      geom_col() +
      coord_flip() +
      labs(
        title = "Pathways Affected Across Multiple Cell Types",
        subtitle = "Female-protective pathways lost or attenuated in T2D",
        x = "Pathway",
        y = "Number of Cell Types",
        fill = "Category"
      ) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 9))
    
    ggsave(file.path(go_dir, "shared_pathways_across_celltypes.png"), 
           p_shared, width = 12, height = 8)
  }
}

cat(sprintf("\nGO analysis complete! Results saved to: %s\n", go_dir))

# ============================================================================
# 8. KEY FINDINGS FOR ABSTRACT
# ============================================================================

cat("\n\n=== KEY FINDINGS FOR ABSTRACT ===\n\n")

# Overall statistics
overall_stats <- combined %>%
  summarise(
    total_genes = n(),
    genes_female_prot_LC = sum(sig_LC & logFC_sex_LC < 0, na.rm = TRUE),
    genes_female_prot_T2D = sum(sig_T2D & logFC_sex_T2D < 0, na.rm = TRUE),
    genes_lost_protection = sum(pattern == "Female-protective lost in T2D", na.rm = TRUE),
    genes_attenuated = sum(protection_attenuated, na.rm = TRUE),
    pct_reduction = 100 * (1 - genes_female_prot_T2D / genes_female_prot_LC)
  )

cat("OVERALL FINDINGS:\n")
cat(sprintf("- Female-protective genes in Lean Controls: %d\n", overall_stats$genes_female_prot_LC))
cat(sprintf("- Female-protective genes in T2D: %d (%.1f%% reduction)\n", 
            overall_stats$genes_female_prot_T2D, overall_stats$pct_reduction))
cat(sprintf("- Genes completely losing female protection: %d\n", overall_stats$genes_lost_protection))
cat(sprintf("- Genes with attenuated protection: %d\n\n", overall_stats$genes_attenuated))

# Cell type specific insights
cat("CELL TYPE-SPECIFIC INSIGHTS:\n")
print(summary_by_celltype %>% dplyr::select(celltype, n_female_prot_LC, n_female_prot_T2D, pct_lost, n_attenuated))

# Top genes with lost/attenuated protection
cat("\nTOP GENES WITH LOST/ATTENUATED FEMALE PROTECTION:\n")
top_genes <- combined %>%
  filter(pattern %in% c("Female-protective lost in T2D", "Attenuated female protection")) %>%
  arrange(desc(abs(delta_logFC))) %>%
  dplyr::select(summary.gene, celltype, logFC_sex_LC, logFC_sex_T2D, delta_logFC, pattern) %>%
  head(20)

print(top_genes)

# Generate GO pathway summary for abstract
cat("\n\nKEY PATHWAY FINDINGS:\n")

# Count pathways by category across all cell types
pathway_summary <- tibble(
  celltype = names(all_go_results),
  n_lc_female_pathways = map_int(all_go_results, ~ifelse(is.null(.x$lc_female), 0, nrow(as.data.frame(.x$lc_female)))),
  n_t2d_female_pathways = map_int(all_go_results, ~ifelse(is.null(.x$t2d_female), 0, nrow(as.data.frame(.x$t2d_female)))),
  n_lost_pathways = map_int(all_go_results, ~ifelse(is.null(.x$lost), 0, nrow(as.data.frame(.x$lost)))),
  n_attenuated_pathways = map_int(all_go_results, ~ifelse(is.null(.x$attenuated), 0, nrow(as.data.frame(.x$attenuated))))
) %>%
  mutate(
    pct_pathway_reduction = 100 * (1 - n_t2d_female_pathways / n_lc_female_pathways)
  )

print(pathway_summary)
write_csv(pathway_summary, file.path(go_dir, "pathway_summary_by_celltype.csv"))

# Highlight interesting shared pathways
if (exists("pathway_counts") && nrow(pathway_counts) > 0) {
  cat("\nPathways appearing in multiple cell types:\n")
  print(pathway_counts %>% head(10))
}

# Save all results
write_csv(combined, file.path(results_dir, "sex_comparison_all_genes.csv"))
write_csv(summary_by_celltype, file.path(results_dir, "sex_comparison_summary_by_celltype.csv"))
write_csv(top_genes, file.path(results_dir, "top_genes_attenuated_protection.csv"))

cat("\n=== ANALYSIS COMPLETE ===\n")
cat(sprintf("Results saved to: %s\n", results_dir))
cat("  - CSV files with all results\n")
cat("  - PNG plots for each cell type\n")








###GSEA 

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)

# ============================================================================
# SETUP
# ============================================================================

# Define directories
base_dir <- "C:/Users/netio/Documents/UofW/Projects/Sex_based_Analysis"
lc_dir <- file.path(base_dir, "LeanControl_Only")
t2d_dir <- file.path(base_dir, "T2D_Only")
gsea_results_dir <- file.path(base_dir, "GSEA_Results")

# Create results directory
if (!dir.exists(gsea_results_dir)) {
  dir.create(gsea_results_dir)
}

cat(sprintf("Results will be saved to: %s\n\n", gsea_results_dir))

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading data...\n")

# Get all files
lc_files <- list.files(lc_dir, pattern = "Full_NEBULA_.*_LC_pooledoffset\\.csv", full.names = TRUE)
t2d_files <- list.files(t2d_dir, pattern = "Full_NEBULA_.*_T2D_pooledoffset\\.csv", full.names = TRUE)

# Extract cell type from filename
extract_celltype <- function(filepath) {
  filename <- basename(filepath)
  gsub("Full_NEBULA_(.+?)_(LC|T2D)_pooledoffset\\.csv", "\\1", filename)
}

# Load and combine data
load_data <- function(files, condition) {
  map_dfr(files, function(f) {
    read_csv(f, show_col_types = FALSE) %>%
      mutate(
        celltype = extract_celltype(f),
        condition = condition
      )
  })
}

lc_data <- load_data(lc_files, "LC")
t2d_data <- load_data(t2d_files, "T2D")

# Prepare sex effect data
prepare_sex_data <- function(df) {
  df %>%
    dplyr::select(
      gene = summary.gene,
      celltype,
      condition,
      logFC = summary.logFC_sexMale,
      pval = summary.p_sexMale
    ) %>%
    # Remove duplicates if any
    distinct(gene, celltype, condition, .keep_all = TRUE) %>%
    # Remove genes with missing values
    filter(!is.na(logFC), !is.na(pval), is.finite(logFC))
}

lc_sex <- prepare_sex_data(lc_data)
t2d_sex <- prepare_sex_data(t2d_data)

# Get all unique cell types
celltypes <- intersect(unique(lc_sex$celltype), unique(t2d_sex$celltype))
cat(sprintf("Found %d cell types: %s\n\n", length(celltypes), paste(celltypes, collapse = ", ")))

# ============================================================================
# GSEA FUNCTION
# ============================================================================

run_gsea <- function(gene_list, ont_type = "BP", analysis_name = "") {
  
  # Remove any NA or infinite values
  gene_list <- gene_list[is.finite(gene_list)]
  
  if (length(gene_list) < 100) {
    cat(sprintf("  Skipping %s - too few genes (%d)\n", analysis_name, length(gene_list)))
    return(NULL)
  }
  
  tryCatch({
    gsea_result <- gseGO(
      geneList = gene_list,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = ont_type,  # BP, CC, or MF
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      by = "fgsea"
    )
    
    if (!is.null(gsea_result) && nrow(as.data.frame(gsea_result)) > 0) {
      cat(sprintf("  ✓ %s (%s): %d pathways\n", analysis_name, ont_type, nrow(as.data.frame(gsea_result))))
      return(gsea_result)
    } else {
      cat(sprintf("  - %s (%s): No significant pathways\n", analysis_name, ont_type))
      return(NULL)
    }
  }, error = function(e) {
    cat(sprintf("  ✗ %s (%s): Error - %s\n", analysis_name, ont_type, e$message))
    return(NULL)
  })
}

# ============================================================================
# ANALYZE EACH CELL TYPE
# ============================================================================

all_results <- list()
all_comparisons <- list()

for (ct in celltypes) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("ANALYZING: %s\n", ct))
  cat(sprintf("========================================\n"))
  
  # Create cell type directory
  ct_dir <- file.path(gsea_results_dir, ct)
  if (!dir.exists(ct_dir)) {
    dir.create(ct_dir)
  }
  
  # Get data for this cell type
  lc_ct <- lc_sex %>% filter(celltype == ct)
  t2d_ct <- t2d_sex %>% filter(celltype == ct)
  
  # Create ranked gene lists (rank by -log10(pval) * sign(logFC))
  create_ranked_list <- function(df) {
    ranked <- df %>%
      mutate(
        rank_metric = -log10(pval) * sign(logFC)
      ) %>%
      arrange(desc(rank_metric)) %>%
      dplyr::select(gene, rank_metric)
    
    # Create named vector
    gene_list <- ranked$rank_metric
    names(gene_list) <- ranked$gene
    
    # Remove duplicates (keep first)
    gene_list <- gene_list[!duplicated(names(gene_list))]
    
    return(gene_list)
  }
  
  lc_ranked <- create_ranked_list(lc_ct)
  t2d_ranked <- create_ranked_list(t2d_ct)
  
  cat(sprintf("LC ranked genes: %d\n", length(lc_ranked)))
  cat(sprintf("T2D ranked genes: %d\n", length(t2d_ranked)))
  
  # Save ranked gene lists
  write_csv(
    tibble(gene = names(lc_ranked), rank_metric = lc_ranked),
    file.path(ct_dir, "LC_ranked_genes.csv")
  )
  write_csv(
    tibble(gene = names(t2d_ranked), rank_metric = t2d_ranked),
    file.path(ct_dir, "T2D_ranked_genes.csv")
  )
  
  # Calculate difference in differences
  # For genes present in both conditions
  common_genes <- intersect(names(lc_ranked), names(t2d_ranked))
  
  delta_ranked <- tibble(
    gene = common_genes,
    lc_rank = lc_ranked[common_genes],
    t2d_rank = t2d_ranked[common_genes],
    delta_rank = t2d_ranked[common_genes] - lc_ranked[common_genes]
  ) %>%
    arrange(desc(abs(delta_rank)))
  
  write_csv(delta_ranked, file.path(ct_dir, "delta_ranked_genes.csv"))
  
  # Create delta ranked list for GSEA
  delta_gene_list <- delta_ranked$delta_rank
  names(delta_gene_list) <- delta_ranked$gene
  
  cat(sprintf("Delta ranked genes: %d\n", length(delta_gene_list)))
  
  # Initialize results storage for this cell type
  all_results[[ct]] <- list()
  
  # ============================================================================
  # RUN GSEA FOR ALL THREE ONTOLOGIES
  # ============================================================================
  
  for (ont in c("BP", "CC", "MF")) {
    cat(sprintf("\n--- %s Ontology ---\n", ont))
    
    # Run GSEA on LC
    gsea_lc <- run_gsea(
      lc_ranked, ont,
      paste(ct, "LC")
    )
    
    # Run GSEA on T2D
    gsea_t2d <- run_gsea(
      t2d_ranked, ont,
      paste(ct, "T2D")
    )
    
    # Run GSEA on delta (difference in differences)
    gsea_delta <- run_gsea(
      delta_gene_list, ont,
      paste(ct, "Delta (T2D - LC)")
    )
    
    # Store results
    all_results[[ct]][[ont]] <- list(
      lc = gsea_lc,
      t2d = gsea_t2d,
      delta = gsea_delta
    )
    
    # Save CSV files
    ont_dir <- file.path(ct_dir, ont)
    if (!dir.exists(ont_dir)) {
      dir.create(ont_dir)
    }
    
    if (!is.null(gsea_lc)) {
      write_csv(as.data.frame(gsea_lc), 
                file.path(ont_dir, "GSEA_LC.csv"))
    }
    if (!is.null(gsea_t2d)) {
      write_csv(as.data.frame(gsea_t2d), 
                file.path(ont_dir, "GSEA_T2D.csv"))
    }
    if (!is.null(gsea_delta)) {
      write_csv(as.data.frame(gsea_delta), 
                file.path(ont_dir, "GSEA_Delta.csv"))
    }
    
    # ============================================================================
    # CREATE VISUALIZATIONS
    # ============================================================================
    
    # GSEA plots for LC
    if (!is.null(gsea_lc) && nrow(as.data.frame(gsea_lc)) > 0) {
      # Dot plot
      p <- dotplot(gsea_lc, showCategory = 20, 
                   title = paste(ct, "- LC Sex Differences (", ont, ")"))
      ggsave(file.path(ont_dir, "GSEA_LC_dotplot.png"), p, width = 12, height = 8)
      
      # Enrichment plot for top pathways
      if (nrow(as.data.frame(gsea_lc)) >= 5) {
        top_pathways <- as.data.frame(gsea_lc) %>% 
          arrange(pvalue) %>% 
          head(5) %>% 
          pull(ID)
        
        for (i in seq_along(top_pathways)) {
          p <- gseaplot2(gsea_lc, geneSetID = top_pathways[i], 
                         title = paste(ct, "LC -", ont))
          ggsave(file.path(ont_dir, paste0("GSEA_LC_pathway_", i, ".png")), 
                 p, width = 10, height = 6)
        }
      }
    }
    
    # GSEA plots for T2D
    if (!is.null(gsea_t2d) && nrow(as.data.frame(gsea_t2d)) > 0) {
      # Dot plot
      p <- dotplot(gsea_t2d, showCategory = 20,
                   title = paste(ct, "- T2D Sex Differences (", ont, ")"))
      ggsave(file.path(ont_dir, "GSEA_T2D_dotplot.png"), p, width = 12, height = 8)
      
      # Enrichment plot for top pathways
      if (nrow(as.data.frame(gsea_t2d)) >= 5) {
        top_pathways <- as.data.frame(gsea_t2d) %>% 
          arrange(pvalue) %>% 
          head(5) %>% 
          pull(ID)
        
        for (i in seq_along(top_pathways)) {
          p <- gseaplot2(gsea_t2d, geneSetID = top_pathways[i],
                         title = paste(ct, "T2D -", ont))
          ggsave(file.path(ont_dir, paste0("GSEA_T2D_pathway_", i, ".png")), 
                 p, width = 10, height = 6)
        }
      }
    }
    
    # GSEA plots for Delta (difference in differences)
    if (!is.null(gsea_delta) && nrow(as.data.frame(gsea_delta)) > 0) {
      # Dot plot
      p <- dotplot(gsea_delta, showCategory = 20,
                   title = paste(ct, "- Change in Sex Differences (T2D - LC) (", ont, ")"))
      ggsave(file.path(ont_dir, "GSEA_Delta_dotplot.png"), p, width = 12, height = 8)
      
      # Enrichment plot for top pathways
      if (nrow(as.data.frame(gsea_delta)) >= 5) {
        top_pathways <- as.data.frame(gsea_delta) %>% 
          arrange(pvalue) %>% 
          head(5) %>% 
          pull(ID)
        
        for (i in seq_along(top_pathways)) {
          p <- gseaplot2(gsea_delta, geneSetID = top_pathways[i],
                         title = paste(ct, "Delta -", ont))
          ggsave(file.path(ont_dir, paste0("GSEA_Delta_pathway_", i, ".png")), 
                 p, width = 10, height = 6)
        }
      }
    }
    
    # ============================================================================
    # COMPARE LC VS T2D PATHWAYS
    # ============================================================================
    
    if (!is.null(gsea_lc) && !is.null(gsea_t2d)) {
      lc_df <- as.data.frame(gsea_lc) %>%
        dplyr::select(ID, Description, NES, pvalue, p.adjust) %>%
        mutate(condition = "LC")
      
      t2d_df <- as.data.frame(gsea_t2d) %>%
        dplyr::select(ID, Description, NES, pvalue, p.adjust) %>%
        mutate(condition = "T2D")
      
      # Merge and compare NES scores
      comparison <- full_join(
        lc_df %>% dplyr::select(ID, Description, NES_LC = NES, padj_LC = p.adjust),
        t2d_df %>% dplyr::select(ID, Description, NES_T2D = NES, padj_T2D = p.adjust),
        by = c("ID", "Description")
      ) %>%
        mutate(
          in_LC = !is.na(NES_LC),
          in_T2D = !is.na(NES_T2D),
          NES_LC = replace_na(NES_LC, 0),
          NES_T2D = replace_na(NES_T2D, 0),
          delta_NES = NES_T2D - NES_LC,
          pattern = case_when(
            in_LC & in_T2D ~ "Shared",
            in_LC & !in_T2D ~ "LC only",
            !in_LC & in_T2D ~ "T2D only",
            TRUE ~ "Neither"
          )
        ) %>%
        arrange(desc(abs(delta_NES)))
      
      all_comparisons[[paste(ct, ont, sep = "_")]] <- comparison
      
      write_csv(comparison, file.path(ont_dir, "LC_vs_T2D_NES_comparison.csv"))
      
      # Visualization 1: Scatter plot of NES scores
      if (nrow(comparison %>% filter(in_LC | in_T2D)) > 0) {
        p <- comparison %>%
          filter(in_LC | in_T2D) %>%
          ggplot(aes(x = NES_LC, y = NES_T2D, color = pattern)) +
          geom_point(alpha = 0.6, size = 3) +
          geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
          geom_hline(yintercept = 0, linetype = "dotted") +
          geom_vline(xintercept = 0, linetype = "dotted") +
          labs(
            title = paste(ct, "- NES Comparison (", ont, ")"),
            subtitle = "LC vs T2D",
            x = "Normalized Enrichment Score (LC)",
            y = "Normalized Enrichment Score (T2D)",
            color = "Pattern",
            caption = "Points above diagonal = stronger enrichment in T2D\nPoints below diagonal = stronger enrichment in LC"
          ) +
          theme_minimal() +
          theme(legend.position = "bottom")
        
        ggsave(file.path(ont_dir, "LC_vs_T2D_NES_scatter.png"), p, width = 10, height = 8)
      }
      
      # Visualization 2: Top pathways with largest change
      if (nrow(comparison %>% filter(in_LC & in_T2D)) > 0) {
        p <- comparison %>%
          filter(in_LC & in_T2D) %>%
          arrange(desc(abs(delta_NES))) %>%
          head(20) %>%
          mutate(Description = str_wrap(Description, 50)) %>%
          ggplot(aes(x = reorder(Description, abs(delta_NES)), y = delta_NES)) +
          geom_col(aes(fill = delta_NES > 0)) +
          coord_flip() +
          scale_fill_manual(values = c("steelblue", "coral"), 
                            labels = c("Stronger in LC", "Stronger in T2D")) +
          labs(
            title = paste(ct, "- Pathways with Largest Change (", ont, ")"),
            subtitle = "Change in Normalized Enrichment Score (T2D - LC)",
            x = "Pathway",
            y = "Delta NES (T2D - LC)",
            fill = ""
          ) +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 9))
        
        ggsave(file.path(ont_dir, "LC_vs_T2D_delta_NES.png"), p, width = 12, height = 10)
      }
    }
  }
}

# ============================================================================
# CROSS-CELL TYPE SUMMARY
# ============================================================================

cat("\n\n========================================\n")
cat("CREATING CROSS-CELL TYPE SUMMARY\n")
cat("========================================\n")

# Create summary table
summary_table <- map_dfr(celltypes, function(ct) {
  bp_lc <- all_results[[ct]][["BP"]]$lc
  bp_t2d <- all_results[[ct]][["BP"]]$t2d
  bp_delta <- all_results[[ct]][["BP"]]$delta
  
  tibble(
    celltype = ct,
    BP_LC_pathways = ifelse(is.null(bp_lc), 0, nrow(as.data.frame(bp_lc))),
    BP_T2D_pathways = ifelse(is.null(bp_t2d), 0, nrow(as.data.frame(bp_t2d))),
    BP_Delta_pathways = ifelse(is.null(bp_delta), 0, nrow(as.data.frame(bp_delta)))
  )
})

write_csv(summary_table, file.path(gsea_results_dir, "summary_by_celltype.csv"))
print(summary_table)

# Find pathways with consistent changes across cell types
for (ont in c("BP", "CC", "MF")) {
  cat(sprintf("\n--- %s Ontology Cross-Cell Type Analysis ---\n", ont))
  
  # Collect all delta GSEA results
  delta_pathways <- map_dfr(celltypes, function(ct) {
    gsea_delta <- all_results[[ct]][[ont]]$delta
    if (!is.null(gsea_delta)) {
      as.data.frame(gsea_delta) %>%
        mutate(celltype = ct) %>%
        dplyr::select(celltype, Description, NES, pvalue, p.adjust)
    }
  })
  
  if (nrow(delta_pathways) > 0) {
    # Find pathways appearing in multiple cell types
    pathway_freq <- delta_pathways %>%
      group_by(Description) %>%
      summarise(
        n_celltypes = n(),
        mean_NES = mean(NES),
        consistent_direction = all(sign(NES) == sign(NES[1]))
      ) %>%
      filter(n_celltypes > 1) %>%
      arrange(desc(n_celltypes), desc(abs(mean_NES)))
    
    if (nrow(pathway_freq) > 0) {
      cat(sprintf("  Found %d pathways in multiple cell types\n", nrow(pathway_freq)))
      
      write_csv(pathway_freq, 
                file.path(gsea_results_dir, paste0("shared_delta_pathways_", ont, ".csv")))
      
      # Visualize
      p <- pathway_freq %>%
        head(20) %>%
        ggplot(aes(x = reorder(str_wrap(Description, 40), n_celltypes), 
                   y = n_celltypes,
                   fill = mean_NES)) +
        geom_col() +
        coord_flip() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
        labs(
          title = paste("Pathways with Changed Sex Differences Across Cell Types (", ont, ")"),
          subtitle = "Based on Delta GSEA (T2D - LC)",
          x = "Pathway",
          y = "Number of Cell Types",
          fill = "Mean NES"
        ) +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 9))
      
      ggsave(file.path(gsea_results_dir, paste0("shared_delta_pathways_", ont, ".png")), 
             p, width = 12, height = 8)
    }
  }
}

cat("\n\n========================================\n")
cat("GSEA ANALYSIS COMPLETE!\n")
cat(sprintf("Results saved to: %s\n", gsea_results_dir))
cat("========================================\n")
