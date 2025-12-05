library(arsenal)
library(Biobase)
library(BiocGenerics)
library(BiocParallel)
library(broom.mixed)
library(colorspace)
library(cowplot)
library(data.table)
library(DirichletReg)
library(dplyr)
library(edgeR)
library(emmeans)
library(enrichR)
library(foreach)
library(future)
library(future.apply)
library(GSEABase)
library(ggdendro)
library(ggpubr)
library(glmmTMB)
library(harmony)
library(jsonlite)
library(kableExtra)
library(limma)
library(MAST)
library(Matrix)
library(msigdbr)
library(muscat)
library(NMF)
library(nebula)
library(patchwork)
library(pheatmap)
library(readxl)
library(REDCapR)
library(reshape2)
library(rstatix)
library(SAVER)
# library(scDC)
library(scater)
library(scran)
library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(tidyverse)
library(UpSetR)
library(WriteXLS)
library(parallel)
library(doParallel)
library(quantreg)
library(aws.s3)

# Set up environment for Kopah
keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

# Pull in PB90 Seurat object
pb90 <- s3readRDS(object = "Kidney transcriptomics/Single cell RNA seq/PB_90_RPCAFix_Metadata_LR_RM_kpmpV1labelled.rds", 
                  bucket = "scrna",  region = "")

# Pull in clinical dataset
rh_rh2_croc_improve_unique <- 
  s3readRDS(object = "data_clean/rh_rh2_croc_improve_unique.RDS", 
            bucket = "harmonized.dataset",
            region = "")

# subset to T2D and OB from RH/RH2/IMPROVE Baseline, HC from CROC
pb90_subset <- subset(pb90, group != "Type_1_Diabetes" & T2D_HC_Phil != "HC_igA" & visit == "baseline")
length(unique(pb90_subset$record_id)) # 49 people with biopsies

# add metadata to PB90 subset
meta_subset <- left_join(subset(pb90_subset@meta.data, select = -c(age, eGFR_CKD_epi, hba1c)), 
                         rh_rh2_croc_improve_unique,
                         by = "kit_id") 
length(unique(meta_subset$record_id.x)) # matches 49 people with biopsies
rownames(meta_subset) <- meta_subset$barcode
pb90_subset@meta.data <- meta_subset

# add general cell type labels
celltype_groups <- list(
  PT = c("PT-S1/S2", "PT-S3", "aPT"),
  TAL = c("C-TAL-1", "C-TAL-2", "aTAL", "dTAL"),
  PC = c("CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC"),
  IC = c("IC-A", "IC-B", "aIC"),
  DTL_ATL = c("DTL", "aDTL", "ATL"),   # grouped thin limbs
  DCT_CNT = c("DCT", "dDCT", "CNT"),   # grouped distal tubule/connecting tubule
  EC = c("EC-AVR", "EC-GC", "EC-PTC", "EC-AEA", "EC-LYM", "EC/VSMC"), 
  Immune = c("MAC", "MON", "cDC", "pDC", "CD4+ T", "CD8+ T", "B", "NK", "cycT"),
  # Immune_Myeloid = c("MAC", "MON", "cDC", "pDC"),
  # Immune_Lymphoid = c("CD4+ T", "CD8+ T", "B", "NK", "cycT"),
  VSMC_P_FIB = c("VSMC/P", "FIB"),
  POD = "POD",
  MC = "MC",                         # mesangial cells
  PEC = "PEC",                       # parietal epithelial cells
  Schwann = "SchwannCells",
  Other = c("non-specific")          # catchall
)

map_celltype_to_general <- function(celltype, celltype_groups) {
  for (group_name in names(celltype_groups)) {
    if (celltype %in% celltype_groups[[group_name]]) {
      return(group_name)
    }
  }
  return("Other")  # For any celltype not found in groups
}

pb90_subset$KPMP_celltype_general <- sapply(pb90_subset$KPMP_celltype, 
                                            map_celltype_to_general, 
                                            celltype_groups = celltype_groups)

# filter out rare gene
ncol(pb90_subset) # Number of cells before filtering (117859)
nrow(pb90_subset) # Number of genes before filtering (31332)

# Filter out rare genes expressed in less than "gene_pct" of cells
expr_matrix <- as.matrix(GetAssayData(pb90_subset, assay = "RNA", layer = "counts"))

# Get cell type information (adjust column name if needed)
cell_types <- pb90_subset@meta.data$KPMP_celltype_general

# Calculate proportion of cells expressing each gene within each cell type
gene_pct_threshold <- 0.05  # 5% threshold

# Initialize a matrix to store proportions for each gene in each cell type
unique_cell_types <- unique(cell_types)
gene_proportions_by_type <- matrix(0, 
                                   nrow = nrow(expr_matrix), 
                                   ncol = length(unique_cell_types))
rownames(gene_proportions_by_type) <- rownames(expr_matrix)
colnames(gene_proportions_by_type) <- unique_cell_types

# Calculate proportion for each cell type
for(ct in unique_cell_types) {
  # Get cells of this type
  cells_of_type <- which(cell_types == ct)
  
  # Calculate proportion of cells expressing each gene in this cell type
  expr_subset <- expr_matrix[, cells_of_type]
  num_cells_expressing <- rowSums(expr_subset > 0)
  total_cells_in_type <- length(cells_of_type)
  
  gene_proportions_by_type[, ct] <- num_cells_expressing / total_cells_in_type
}

# Keep genes expressed in at least 5% of cells in AT LEAST ONE cell type
max_proportion_per_gene <- apply(gene_proportions_by_type, 1, max)
genes_to_keep <- names(max_proportion_per_gene[max_proportion_per_gene >= gene_pct_threshold])

# Check results
cat("Total genes:", nrow(expr_matrix), "\n") # Total genes: 31332 
cat("Genes to keep:", length(genes_to_keep), "\n") # Genes to keep: 16286 
cat("Proportion kept:", round(length(genes_to_keep)/nrow(expr_matrix)*100, 1), "%\n") # Proportion kept: 52 %

# Optional: Check distribution of max proportions
hist(max_proportion_per_gene)

# Filter the Seurat object
pb90_subset <- subset(pb90_subset, features = genes_to_keep)


#Check the number of Mitochondrial genes to start
sum(grepl("^MT-", rownames(pb90_subset))) 

# Identify mitochondrial genes (human: start with "MT-")
mito_genes <- grep("^MT-", rownames(pb90_subset), value = TRUE)
pb90_subset <- subset(pb90_subset, features = setdiff(rownames(pb90_subset), mito_genes))

#Check the number of Mitochondrial genes after filtering to ensure filtering step was successful
sum(grepl("^MT-", rownames(pb90_subset))) #Should be 0

# Identify ribosomal genes
ribo_genes <- c(
  "RPL22", "RPL11", "RPS8", "RPL5", "RPS27", "RPS7", "RPS27A", "RPL31", "RPL37A", "RPL32", "RPL15", "RPL14", "RPL29",
  "RPL24", "RPL22L1", "RPL35A", "RPL9", "RPL34", "RPS3A", "RPL37", "RPS23", "RPS14", "RPS18", "RPS10", "RPL10A", 
  "RPS20", "RPL7", "RPL30", "RPL8", "RPS6", "RPL35", "RPL12", "RPL7A", "RPS24", "RPLP2", "RPL27A", "RPS13", "RPS3",
  "RPS25", "RPS26", "RPL41", "RPL6", "RPLP0", "RPL21", "RPS29", "RPL4", "RPLP1", "RPS17", "RPS2", "RPS15A", "RPL13",
  "RPL26", "RPL23A", "RPL23", "RPL19", "RPL27", "RPL38", "RPL17", "RPS15", "RPL36", "RPS28", "RPL18A", "RPS16", 
  "RPS19", "RPL18", "RPL13A", "RPS11", "RPS9", "RPL28", "RPS5", "RPS21", "RPL3", "RPS4X", "RPL36A", "RPL39", 
  "RPL10", "RPS4Y1"
) # grep("^RPL|^RPS", rownames(attempt_so_raw), value = TRUE) captures some none ribosomal genes

pb90_subset <- subset(pb90_subset, features = setdiff(rownames(pb90_subset), ribo_genes))


# Check if data are normalized
Assays(pb90_subset)  # Should show e.g., "RNA" with "counts" and "data"
pb90_subset@assays # Another method to see RNA assays
head(GetAssayData(pb90_subset, layer = "counts")[, 1:5])  # Raw counts
head(GetAssayData(pb90_subset, layer = "data")[, 1:5])    # Normalized data


# Exploration of data
# Randomly select 100 genes from the Seurat object
genes <- sample(rownames(pb90_subset), 100)

# Set up a 2x3 plotting layout so you can plot multiple histograms in one figure
par(mfrow = c(2, 3))

# Loop over each randomly selected gene
for (g in genes) {
  
  # Plot a histogram of the expression values for gene 'g'
  # using normalized expression values from the "data" slot
  hist(GetAssayData(pb90_subset, layer = "data")[g, ],
       main = g,                # Title of the plot = gene name
       xlab = "Normalized Expression")     # Label for x-axis
  
  # using raw counts expression values from the "counts" slot
  hist(GetAssayData(pb90_subset, layer = "counts")[g, ],
       main = g,                # Title of the plot = gene name
       xlab = "Raw Counts Expression")     # Label for x-axis
}

# Save Seurat object for downstream analysis
s3saveRDS(pb90_subset, object = "data_clean/subset/pb90_ckd_analysis_subset.rds", bucket = "scrna", region = "")

# Subset Seurat object into general cell types and save
for (cell in names(celltype_groups)) {
  gc()
  pb90_celltype <- subset(pb90_subset, 
                          KPMP_celltype_general == cell)
  s3saveRDS(pb90_celltype, object = paste0("data_clean/subset/pb90_ckd_analysis/pb90_ckd_analysis_subset_", cell, ".rds"), 
            bucket = "scrna", region = "", multipart = T)
  
}

pb90_subset <- s3readRDS(object = "data_clean/subset/pb90_ckd_analysis_subset.rds", bucket = "scrna", region = "")

length(unique(pb90_subset$kit_id))

library(ggrepel)
library(RColorBrewer)
library(ggtrace)
my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(43)
# Read Seurat object and run UMAP, obtain cell numbers, etc
pb90_subset_meta <- pb90_subset@meta.data %>%
  cbind(pb90_subset@reductions$umap.harmony@cell.embeddings)

s3saveRDS(pb90_subset_meta, object = "data_clean/subset/pb90_ckd_analysis_subset_meta.rds", bucket = "scrna", region = "")
pb90_subset_meta <- s3readRDS(object = "data_clean/subset/pb90_ckd_analysis_subset_meta.rds", bucket = "scrna", region = "")

centers <- pb90_subset_meta %>%
  group_by(KPMP_celltype) %>%
  summarise(x = median(umapharmony_1),
            y = median(umapharmony_2), .groups = "drop")

umap_p <- ggplot(pb90_subset_meta,
       aes(x = umapharmony_1, y = umapharmony_2, color = KPMP_celltype)) +
  geom_point(alpha = 0.2) +
  geom_text_repel(data = centers,
                  aes(x = x, y = y, label = KPMP_celltype),
                  inherit.aes = FALSE,
                  size = 4, fontface = "bold", color = "black",
                  family = "arial") +
  scale_color_manual(values = my_colors) +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(hjust = 0),   # move x label to left
        axis.title.y = element_text(hjust = 0))  + # push y label down to bottom) 
  labs(x = "UMAP1", y = "UMAP2")

# Save to a temporary connection
tmp <- tempfile(fileext = ".png")
ggsave(tmp, plot = umap_p, width = 10, height = 8, dpi = 300)

# Upload directly
put_object(
  file = tmp,
  object = "Projects/CKD/RH_RH2/Results/Figures/rhrh2imp_umap_kpmpcelltype.png",
  bucket = "scrna",
  region = ""
)

unlink(tmp)

cell_counts <- pb90_subset_meta %>%
  dplyr::count(KPMP_celltype_general) %>%
  arrange(desc(n))

# Create bar plot
cell_count_bar <- ggplot(cell_counts, aes(x = reorder(KPMP_celltype_general, n), y = n)) +
  geom_bar(stat = "identity", fill = "#97a97c") +
  geom_text(aes(label = n), hjust = -0.1, size = 4) +
  coord_flip() +
  labs(x = NULL, 
       y = "Count") +
  theme(text = element_text(size = 15),
        legend.position = "none",
        panel.grid = element_blank()) +
  theme_transparent +
  scale_y_continuous(expand = expansion(add = 5000)) 

# Save to a temporary connection
tmp <- tempfile(fileext = ".png")
ggsave(tmp, plot = cell_count_bar, width = 10, height = 8, dpi = 300)

# Upload directly
put_object(
  file = tmp,
  object = "Projects/CKD/RH_RH2/Results/Figures/rhrh2imp_kpmp_cell_count_bar.png",
  bucket = "scrna",
  region = ""
)

pb90_subset_meta %>%
  filter(KPMP_celltype_general == "PT") %>%
  filter(group != "Lean_Control") %>%
  ggplot(aes(x = record_id, fill = KPMP_celltype)) +
  geom_bar(position = "fill") +
  ylab("Proportion") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw() +
  facet_wrap(~dkd_group_30) +
  theme(axis.text.x = element_text(angle = 60))

pb90_subset_meta %>%
  filter(celltype_harmony %in% c("TAL-1", "TAL-2", "TAL-3")) %>%
  ggplot(aes(x = dkd_group_30, fill = celltype_harmony)) +
  geom_bar(position = "fill") 

pb90_subset_meta %>%
  filter(KPMP_celltype_general == "TAL") %>%
  ggplot(aes(x = dkd_group_30, fill = KPMP_celltype)) +
  geom_bar(position = "fill") +
  ylab("Proportion") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw()

pb90_subset_meta %>%
  filter(KPMP_celltype_general == "EC") %>%
  ggplot(aes(x = dkd_group_30, fill = KPMP_celltype)) +
  geom_bar(position = "fill") +
  ylab("Proportion") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw()

pb90_subset_meta %>%
  filter(KPMP_celltype_general == "Immune") %>%
  ggplot(aes(x = dkd_group_30, fill = KPMP_celltype)) +
  geom_bar(position = "fill") +
  ylab("Proportion") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw()

# Cell type proportion tests
library(speckle)

celltype_general_map <- list(
  PT = c("PT-1", "PT-2", "PT-3", "PT-4", "PT-5", "PT_lowQuality"),
  TAL = c("TAL-1", "TAL-2", "TAL-3", "TAL_highUMI"),
  DTL_ATL = c("DTL-1", "DTL-2", "ATL"),
  DCT_CNT = c("DCT", "CNT"),
  PC = c("PC-1", "PC-2", "tPC-IC"),
  IC = c("IC-A", "IC-B", "IC-A_lowQuality"),
  
  Immune = c("T", "NKT/NKC", "B", "MON", "MAC"),
  
  EC = c("EC-GC", "EC-PTC", "EC-AEA", "EC-LYM"),
  VSMC_P_FIB = c("VSMC/MC/FIB"),
  POD = c("POD"),
  PEC = c("PEC"),
  Other = c("TAL_highUMI")  # if needed or any low-quality unassigned
)

celltype_map <- unlist(lapply(names(celltype_general_map), function(x) {
  setNames(rep(x, length(celltype_general_map[[x]])), celltype_general_map[[x]])
}))

pb90_subset_meta <- pb90_subset_meta %>%
  mutate(celltype_general_from_harmony = recode(celltype_harmony, !!!celltype_map),
         celltype_general_from_harmony = ifelse(is.na(celltype_general_from_harmony), "Other", 
                                                celltype_general_from_harmony))
table(pb90_subset_meta$celltype_harmony, pb90_subset_meta$celltype_general_from_harmony)

propeller_results <- lapply(names(celltype_general_map), function(cell) {
  pb90_subset_meta_celltype <- subset(pb90_subset_meta, 
                                      celltype_general_from_harmony == cell &
                                        group != "Lean_Control")
  
  pb90_subset_meta_celltype$celltype_harmony <- as.character(pb90_subset_meta_celltype$celltype_harmony)
  
  transformed_props <- getTransformedProps(
    clusters = pb90_subset_meta_celltype$celltype_harmony,
    sample = paste0(pb90_subset_meta_celltype$record_id, "_", 
                    pb90_subset_meta_celltype$dkd_group_30)
  )
  
  dkd_group_30 <- ifelse(grepl("_non_DKD$", colnames(transformed_props$Counts)),
                         "non_DKD", "DKD")
  
  design <- model.matrix(~ 0 + dkd_group_30)
  colnames(design) <- gsub("dkd_group_30", "", colnames(design))
  
  contrasts <- c(DKD = 1, non_DKD = -1)
  
  propeller_res <- propeller.ttest(prop.list = transformed_props,
                                   design = design,
                                   contrasts = contrasts,
                                   robust = TRUE,
                                   trend = FALSE,
                                   sort = TRUE)
  print(propeller_res)
  return(propeller_res)
})

subject_props <- pb90_subset_meta %>%
  group_by(record_id, dkd_group_30, KPMP_celltype) %>%
  summarize(cell_n = n(), .groups = "drop") %>%
  group_by(record_id, dkd_group_30) %>%
  mutate(total_n = sum(cell_n),
         prop = cell_n / total_n) %>%
  ungroup()

library(ggpubr)

ggplot(subject_props, aes(x = dkd_group_30, y = prop, fill = dkd_group_30)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     size = 6) +
  facet_wrap(~ KPMP_celltype, scales = "free_y") +
  ylab("Proportion per subject") +
  xlab("Group") +
  theme_bw() +
  theme(text = element_text(size = 10),
        strip.text = element_text(size = 10),
        legend.position = "none") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.2)))

# add general cell type labels
pb90_subset <- s3readRDS(object = "data_clean/subset/pb90_ckd_analysis_subset.rds", bucket = "scrna", region = "")

celltype_groups2 <- list(
  PT = c("PT-S1/S2", "PT-S3", "aPT"),
  TAL = c("C-TAL-1", "C-TAL-2", "aTAL", "dTAL"),
  PC = c("CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC"),
  IC = c("IC-A", "IC-B", "aIC"),
  DTL_ATL = c("DTL", "aDTL", "ATL"),   # grouped thin limbs
  DCT_CNT = c("DCT", "dDCT", "CNT"),   # grouped distal tubule/connecting tubule
  EC = c("EC-AVR", "EC-GC", "EC-PTC", "EC-AEA", "EC-LYM", "EC/VSMC"), 
  # Immune = c("MAC", "MON", "cDC", "pDC", "CD4+ T", "CD8+ T", "B", "NK", "cycT"),
  Immune_Myeloid = c("MAC", "MON", "cDC", "pDC"),
  Immune_Lymphoid = c("CD4+ T", "CD8+ T", "B", "NK", "cycT"),
  VSMC_P_FIB = c("VSMC/P", "FIB"),
  POD = "POD",
  MC = "MC",                         # mesangial cells
  PEC = "PEC",                       # parietal epithelial cells
  Schwann = "SchwannCells",
  Other = c("non-specific")          # catchall
)

map_celltype_to_general <- function(celltype, celltype_groups2) {
  for (group_name in names(celltype_groups2)) {
    if (celltype %in% celltype_groups2[[group_name]]) {
      return(group_name)
    }
  }
  return("Other")  # For any celltype not found in groups
}

pb90_subset$KPMP_celltype_general2 <- sapply(pb90_subset$KPMP_celltype, 
                                            map_celltype_to_general, 
                                            celltype_groups2 = celltype_groups2)

s3saveRDS(pb90_subset, object = "data_clean/subset/pb90_ckd_analysis_subset.rds", bucket = "scrna", region = "", multipart = T)


#################################################################################################
# temp scratch
folders <- c(
  
  # --- Base DKD vs nonDKD ---
  # "DKD_vs_nonDKD_30" = "dkd30",
  # "DKD_vs_nonDKD_100" = "dkd100",
  # 
  # # --- DKD vs nonDKD with T2D ---
  # "DKD_vs_nonDKD_30_t2d" = "dkd30_t2d",
  # "DKD_vs_nonDKD_100_t2d" = "dkd100_t2d",
  # 
  # # --- DKD vs HC (30 / 100) ---
  # "DKD_30_vs_HC" = "dkd30_hc",
  # "DKD_100_vs_HC" = "dkd100_hc",
  # 
  # # --- DKD vs HC with T2D ---
  # "DKD_30_t2d_vs_HC" = "dkd30_t2d_hc",
  # "DKD_100_t2d_vs_HC" = "dkd100_t2d_hc",
  # 
  # # --- nonDKD vs HC (30 / 100) ---
  # "nonDKD_30_vs_HC" = "nondkd30_hc",
  # "nonDKD_100_vs_HC" = "nondkd100_hc",
  # 
  # # --- nonDKD vs HC with T2D ---
  # "nonDKD_30_t2d_vs_HC" = "nondkd30_t2d_hc",
  # "nonDKD_100_t2d_vs_HC" = "nondkd100_t2d_hc",
  # 
  # # --- GLP (within DKD) ---
  # "DKD_30_GLP_Y_vs_DKD_30_GLP_N" = "dkd_30_glpy_glpn",
  # "DKD_100_GLP_Y_vs_DKD_100_GLP_N" = "dkd_100_glpy_glpn",
  # 
  # # --- GLP (within nonDKD) ---
  # "nonDKD_30_GLP_Y_vs_nonDKD_30_GLP_N" = "nondkd_30_glpy_glpn",
  # "nonDKD_100_GLP_Y_vs_nonDKD_100_GLP_N" = "nondkd_100_glpy_glpn",
  # 
  # # --- GLP N/Y vs HC (nonDKD) ---
  # "nonDKD_30_GLP_N_vs_HC" = "nondkd30_glpn_hc",
  # "nonDKD_30_GLP_Y_vs_HC" = "nondkd30_glpy_hc",
  # "nonDKD_100_GLP_N_vs_HC" = "nondkd100_glpn_hc",
  # "nonDKD_100_GLP_Y_vs_HC" = "nondkd100_glpy_hc",
  # 
  # # --- GLP N/Y vs HC (DKD) ---
  # "DKD_30_GLP_N_vs_HC" = "dkd30_glpn_hc",
  # "DKD_30_GLP_Y_vs_HC" = "dkd30_glpy_hc",
  # "DKD_100_GLP_N_vs_HC" = "dkd100_glpn_hc",
  # "DKD_100_GLP_Y_vs_HC" = "dkd100_glpy_hc",
  # 
  # # --- SGLT2i (within DKD) ---
  # "DKD_30_SGLT2i_Y_vs_DKD_30_SGLT2i_N" = "dkd30_sglt2iy_sglt2in",
  "DKD_100_SGLT2i_Y_vs_DKD_100_SGLT2i_N" = "dkd100_sglt2iy_sglt2in",
  
  # --- SGLT2i (within nonDKD) ---
  "nonDKD_30_SGLT2i_Y_vs_nonDKD_30_SGLT2i_N" = "nondkd30_sglt2iy_sglt2in",
  "nonDKD_100_SGLT2i_Y_vs_nonDKD_100_SGLT2i_N" = "nondkd100_sglt2iy_sglt2in",
  
  # --- SGLT2i N/Y vs HC (nonDKD) ---
  "nonDKD_30_SGLT2i_N_vs_HC" = "nondkd30_sglt2in_hc",
  "nonDKD_30_SGLT2i_Y_vs_HC" = "nondkd30_sglt2iy_hc",
  "nonDKD_100_SGLT2i_N_vs_HC" = "nondkd100_sglt2in_hc",
  "nonDKD_100_SGLT2i_Y_vs_HC" = "nondkd100_sglt2iy_hc",
  
  # --- SGLT2i N/Y vs HC (DKD) ---
  "DKD_30_SGLT2i_N_vs_HC" = "dkd30_sglt2in_hc",
  "DKD_30_SGLT2i_Y_vs_HC" = "dkd30_sglt2iy_hc",
  "DKD_100_SGLT2i_N_vs_HC" = "dkd100_sglt2in_hc",
  "DKD_100_SGLT2i_Y_vs_HC" = "dkd100_sglt2iy_hc"
)

analysis_config <- list(
  # DKD vs nonDKD comparisons (ACR >= 100)
  DKD_vs_nonDKD_100 = list(
    subset_cond = "group != 'Lean_Control' & !is.na(dkd_group_100)",
    group_var = "dkd_group_100", ref_level = "non_DKD",
    pval_col = "p_dkd_group_100DKD", logfc_col = "logFC_dkd_group_100DKD",
    s3_subdir = "DKD_vs_nonDKD_100", file_suffix = "dkd100"
  ),
  DKD_vs_nonDKD_100_t2d = list(
    subset_cond = "group == 'Type_2_Diabetes' & !is.na(dkd_group_100)",
    group_var = "dkd_group_100", ref_level = "non_DKD",
    pval_col = "p_dkd_group_100DKD", logfc_col = "logFC_dkd_group_100DKD",
    s3_subdir = "DKD_vs_nonDKD_100_t2d", file_suffix = "dkd100_t2d"
  ),
  DKD_vs_nonDKD_100_nosglt2i = list(
    subset_cond = "group == 'Type_2_Diabetes' & epic_sglti2_1 == 'No' & !is.na(dkd_group_100)",
    group_var = "dkd_group_100", ref_level = "non_DKD",
    pval_col = "p_dkd_group_100DKD", logfc_col = "logFC_dkd_group_100DKD",
    s3_subdir = "DKD_vs_nonDKD_100_nosglt2i", file_suffix = "dkd100_nosglt2i"
  ),
  
  # DKD vs nonDKD comparisons (ACR >= 30)
  DKD_vs_nonDKD_30 = list(
    subset_cond = "group != 'Lean_Control' & !is.na(dkd_group_30)",
    group_var = "dkd_group_30", ref_level = "non_DKD",
    pval_col = "p_dkd_group_30DKD", logfc_col = "logFC_dkd_group_30DKD",
    s3_subdir = "DKD_vs_nonDKD_30", file_suffix = "dkd30"
  ),
  DKD_vs_nonDKD_30_t2d = list(
    subset_cond = "group == 'Type_2_Diabetes' & !is.na(dkd_group_30)",
    group_var = "dkd_group_30", ref_level = "non_DKD",
    pval_col = "p_dkd_group_30DKD", logfc_col = "logFC_dkd_group_30DKD",
    s3_subdir = "DKD_vs_nonDKD_30_t2d", file_suffix = "dkd30_t2d"
  ),
  DKD_vs_nonDKD_30_nosglt2i = list(
    subset_cond = "group == 'Type_2_Diabetes' & epic_sglti2_1 == 'No' & !is.na(dkd_group_30)",
    group_var = "dkd_group_30", ref_level = "non_DKD",
    pval_col = "p_dkd_group_30DKD", logfc_col = "logFC_dkd_group_30DKD",
    s3_subdir = "DKD_vs_nonDKD_30_nosglt2i", file_suffix = "dkd30_nosglt2i"
  ),
  
  # DKD vs HC comparisons (ACR >= 100)
  DKD_100_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'non_DKD' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC",
    pval_col = "p_dkd_group_100_hcDKD", logfc_col = "logFC_dkd_group_100_hcDKD",
    s3_subdir = "DKD_100_vs_HC", file_suffix = "dkd100_hc"
  ),
  DKD_100_t2d_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'non_DKD' & group != 'Obese_Control' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC",
    pval_col = "p_dkd_group_100_hcDKD", logfc_col = "logFC_dkd_group_100_hcDKD",
    s3_subdir = "DKD_100_t2d_vs_HC", file_suffix = "dkd100_t2d_hc"
  ),
  nonDKD_100_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC",
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_vs_HC", file_suffix = "nondkd100_hc"
  ),
  nonDKD_100_t2d_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & group != 'Obese_Control' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC",
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_t2d_vs_HC", file_suffix = "nondkd100_t2d_hc"
  ),
  
  # DKD vs HC comparisons (ACR >= 30)
  DKD_30_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'non_DKD' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC",
    pval_col = "p_dkd_group_30_hcDKD", logfc_col = "logFC_dkd_group_30_hcDKD",
    s3_subdir = "DKD_30_vs_HC", file_suffix = "dkd30_hc"
  ),
  DKD_30_t2d_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'non_DKD' & group != 'Obese_Control' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC",
    pval_col = "p_dkd_group_30_hcDKD", logfc_col = "logFC_dkd_group_30_hcDKD",
    s3_subdir = "DKD_30_t2d_vs_HC", file_suffix = "dkd30_t2d_hc"
  ),
  nonDKD_30_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC",
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_vs_HC", file_suffix = "nondkd30_hc"
  ),
  nonDKD_30_t2d_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & group != 'Obese_Control' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC",
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_t2d_vs_HC", file_suffix = "nondkd30_t2d_hc"
  ),
  
  # GLP-1RA comparisons
  GLP_N_vs_HC = list(
    subset_cond = "glp_t2dob != 'GLP_Y' & !is.na(glp_t2dob)",
    group_var = "glp_t2dob", ref_level = "HC",
    pval_col = "p_glp_t2dobGLP_N", logfc_col = "logFC_glp_t2dobGLP_N",
    s3_subdir = "GLP_N_vs_HC", file_suffix = "glpn_hc"
  ),
  GLP_Y_vs_GLP_N = list(
    subset_cond = "group != 'Lean_Control' & !is.na(glp_t2dob)",
    group_var = "glp_t2dob", ref_level = "GLP_N",
    pval_col = "p_glp_t2dobGLP_Y", logfc_col = "logFC_glp_t2dobGLP_Y",
    s3_subdir = "GLP_Y_vs_GLP_N", file_suffix = "glpy_glpn"
  ),
  
  # GLP within DKD/nonDKD subgroups (ACR >= 100)
  DKD_100_GLP_Y_vs_DKD_100_GLP_N = list(
    subset_cond = "dkd_group_100_hc == 'DKD' & !is.na(glp_t2dob)",
    group_var = "glp_t2dob", ref_level = "GLP_N",
    pval_col = "p_glp_t2dobGLP_Y", logfc_col = "logFC_glp_t2dobGLP_Y",
    s3_subdir = "DKD_100_GLP_Y_vs_DKD_100_GLP_N", file_suffix = "dkd_100_glpy_glpn"
  ),
  nonDKD_100_GLP_Y_vs_nonDKD_100_GLP_N = list(
    subset_cond = "dkd_group_100_hc == 'non_DKD' & !is.na(glp_t2dob)",
    group_var = "glp_t2dob", ref_level = "GLP_N",
    pval_col = "p_glp_t2dobGLP_Y", logfc_col = "logFC_glp_t2dobGLP_Y",
    s3_subdir = "nonDKD_100_GLP_Y_vs_nonDKD_100_GLP_N", file_suffix = "nondkd_100_glpy_glpn"
  ),
  nonDKD_100_GLP_N_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & glp_t2dob != 'GLP_Y' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC",
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_GLP_N_vs_HC", file_suffix = "nondkd100_glpn_hc"
  ),
  nonDKD_100_GLP_Y_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & glp_t2dob != 'GLP_N' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC",
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_GLP_Y_vs_HC", file_suffix = "nondkd100_glpy_hc"
  ),
  DKD_100_GLP_Y_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'non_DKD' & glp_t2dob != 'GLP_N' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC",
    pval_col = "p_dkd_group_100_hcDKD", logfc_col = "logFC_dkd_group_100_hcDKD",
    s3_subdir = "DKD_100_GLP_Y_vs_HC", file_suffix = "dkd100_glpy_hc"
  ),
  DKD_100_GLP_N_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'non_DKD' & glp_t2dob != 'GLP_Y' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC",
    pval_col = "p_dkd_group_100_hcDKD", logfc_col = "logFC_dkd_group_100_hcDKD",
    s3_subdir = "DKD_100_GLP_N_vs_HC", file_suffix = "dkd100_glpn_hc"
  ),
  
  # GLP within DKD/nonDKD subgroups (ACR >= 30)
  DKD_30_GLP_Y_vs_DKD_30_GLP_N = list(
    subset_cond = "dkd_group_30_hc == 'DKD' & !is.na(glp_t2dob)",
    group_var = "glp_t2dob", ref_level = "GLP_N",
    pval_col = "p_glp_t2dobGLP_Y", logfc_col = "logFC_glp_t2dobGLP_Y",
    s3_subdir = "DKD_30_GLP_Y_vs_DKD_30_GLP_N", file_suffix = "dkd_30_glpy_glpn"
  ),
  nonDKD_30_GLP_Y_vs_nonDKD_30_GLP_N = list(
    subset_cond = "dkd_group_30_hc == 'non_DKD' & !is.na(glp_t2dob)",
    group_var = "glp_t2dob", ref_level = "GLP_N",
    pval_col = "p_glp_t2dobGLP_Y", logfc_col = "logFC_glp_t2dobGLP_Y",
    s3_subdir = "nonDKD_30_GLP_Y_vs_nonDKD_30_GLP_N", file_suffix = "nondkd_30_glpy_glpn"
  ),
  nonDKD_30_GLP_N_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & glp_t2dob != 'GLP_Y' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC",
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_GLP_N_vs_HC", file_suffix = "nondkd30_glpn_hc"
  ),
  nonDKD_30_GLP_Y_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & glp_t2dob != 'GLP_N' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC",
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_GLP_Y_vs_HC", file_suffix = "nondkd30_glpy_hc"
  ),
  DKD_30_GLP_Y_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'non_DKD' & glp_t2dob != 'GLP_N' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC",
    pval_col = "p_dkd_group_30_hcDKD", logfc_col = "logFC_dkd_group_30_hcDKD",
    s3_subdir = "DKD_30_GLP_Y_vs_HC", file_suffix = "dkd30_glpy_hc"
  ),
  DKD_30_GLP_N_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'non_DKD' & glp_t2dob != 'GLP_Y' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC",
    pval_col = "p_dkd_group_30_hcDKD", logfc_col = "logFC_dkd_group_30_hcDKD",
    s3_subdir = "DKD_30_GLP_N_vs_HC", file_suffix = "dkd30_glpn_hc"
  ),
  
  # SGLT2i comparisons
  SGLT2i_N_vs_HC = list(
    subset_cond = "sglt2i_t2dob != 'SGLT2i_Y' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_N", logfc_col = "logFC_sglt2i_t2dobSGLT2i_N",
    s3_subdir = "SGLT2i_N_vs_HC", file_suffix = "sglt2in_hc"
  ),
  SGLT2i_Y_vs_SGLT2i_N = list(
    subset_cond = "group != 'Lean_Control' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "SGLT2i_N", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_Y", logfc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    s3_subdir = "SGLT2i_Y_vs_SGLT2i_N", file_suffix = "sglt2iy_sglt2in"
  ),
  DKD_30_SGLT2i_N_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'non_DKD' & sglt2i_t2dob != 'SGLT2i_Y' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_N", logfc_col = "logFC_sglt2i_t2dobSGLT2i_N",
    s3_subdir = "DKD_30_SGLT2i_N_vs_HC", file_suffix = "dkd30_sglt2in_hc"
  ),
  DKD_30_SGLT2i_Y_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'non_DKD' & sglt2i_t2dob != 'SGLT2i_N' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_Y", logfc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    s3_subdir = "DKD_30_SGLT2i_Y_vs_HC", file_suffix = "dkd30_sglt2iy_hc"
  ),
  DKD_30_SGLT2i_Y_vs_DKD_30_SGLT2i_N = list(
    subset_cond = "dkd_group_30_hc == 'DKD' & group != 'Lean_Control' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "SGLT2i_N", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_Y", logfc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    s3_subdir = "DKD_30_SGLT2i_Y_vs_DKD_30_SGLT2i_N", file_suffix = "dkd30_sglt2iy_sglt2in"
  ),
  DKD_100_SGLT2i_Y_vs_DKD_100_SGLT2i_N = list(
    subset_cond = "dkd_group_100_hc == 'DKD' & group != 'Lean_Control' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "SGLT2i_N", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_Y", logfc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    s3_subdir = "DKD_100_SGLT2i_Y_vs_DKD_100_SGLT2i_N", file_suffix = "dkd100_sglt2iy_sglt2in"
  ),
  nonDKD_100_SGLT2i_N_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & sglt2i_t2dob != 'SGLT2i_Y' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_SGLT2i_N_vs_HC", file_suffix = "nondkd100_sglt2in_hc"
  ),
  nonDKD_100_SGLT2i_Y_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & sglt2i_t2dob != 'SGLT2i_N' & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_SGLT2i_Y_vs_HC", file_suffix = "nondkd100_sglt2iy_hc"
  ),
  nonDKD_30_SGLT2i_N_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & sglt2i_t2dob != 'SGLT2i_Y' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_SGLT2i_N_vs_HC", file_suffix = "nondkd30_sglt2in_hc"
  ),
  nonDKD_30_SGLT2i_Y_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & sglt2i_t2dob != 'SGLT2i_N' & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_SGLT2i_Y_vs_HC", file_suffix = "nondkd30_sglt2iy_hc"
  ),
  nonDKD_30_SGLT2i_Y_vs_nonDKD_30_SGLT2i_N = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & group != 'Lean_Control' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "SGLT2i_N", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_Y", logfc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    s3_subdir = "nonDKD_30_SGLT2i_Y_vs_nonDKD_30_SGLT2i_N", file_suffix = "nondkd30_sglt2iy_sglt2in"
  ),
  nonDKD_100_SGLT2i_Y_vs_nonDKD_100_SGLT2i_N = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & group != 'Lean_Control' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "SGLT2i_N", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_Y", logfc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    s3_subdir = "nonDKD_100_SGLT2i_Y_vs_nonDKD_100_SGLT2i_N", file_suffix = "nondkd100_sglt2iy_sglt2in"
  ),
  DKD_100_SGLT2i_N_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'non_DKD' & sglt2i_t2dob != 'SGLT2i_Y' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_N", logfc_col = "logFC_sglt2i_t2dobSGLT2i_N",
    s3_subdir = "DKD_100_SGLT2i_N_vs_HC", file_suffix = "dkd100_sglt2in_hc"
  ),
  DKD_100_SGLT2i_Y_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'non_DKD' & sglt2i_t2dob != 'SGLT2i_N' & !is.na(sglt2i_t2dob)",
    group_var = "sglt2i_t2dob", ref_level = "HC", needs_sglt2i = TRUE,
    pval_col = "p_sglt2i_t2dobSGLT2i_Y", logfc_col = "logFC_sglt2i_t2dobSGLT2i_Y",
    s3_subdir = "DKD_100_SGLT2i_Y_vs_HC", file_suffix = "dkd100_sglt2iy_hc"
  ),
  
  # Combo comparisons - each treatment group vs HC
  Neither_vs_HC = list(
    subset_cond = "combo_t2dob %in% c('Neither', 'HC') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_combo_t2dobNeither", logfc_col = "logFC_combo_t2dobNeither",
    s3_subdir = "Neither_vs_HC", file_suffix = "neither_hc"
  ),
  SGLT2i_only_vs_HC = list(
    subset_cond = "combo_t2dob %in% c('SGLT2i_only', 'HC') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_combo_t2dobSGLT2i_only", logfc_col = "logFC_combo_t2dobSGLT2i_only",
    s3_subdir = "SGLT2i_only_vs_HC", file_suffix = "sglt2ionly_hc"
  ),
  GLP_only_vs_HC = list(
    subset_cond = "combo_t2dob %in% c('GLP_only', 'HC') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_combo_t2dobGLP_only", logfc_col = "logFC_combo_t2dobGLP_only",
    s3_subdir = "GLP_only_vs_HC", file_suffix = "glponly_hc"
  ),
  Both_vs_HC = list(
    subset_cond = "combo_t2dob %in% c('Both', 'HC') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_combo_t2dobBoth", logfc_col = "logFC_combo_t2dobBoth",
    s3_subdir = "Both_vs_HC", file_suffix = "both_hc"
  ),
  
  # Combo comparisons - treatment groups vs Neither (no treatment)
  SGLT2i_only_vs_Neither = list(
    subset_cond = "combo_t2dob %in% c('SGLT2i_only', 'Neither') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "Neither", needs_combo = TRUE,
    pval_col = "p_combo_t2dobSGLT2i_only", logfc_col = "logFC_combo_t2dobSGLT2i_only",
    s3_subdir = "SGLT2i_only_vs_Neither", file_suffix = "sglt2ionly_neither"
  ),
  GLP_only_vs_Neither = list(
    subset_cond = "combo_t2dob %in% c('GLP_only', 'Neither') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "Neither", needs_combo = TRUE,
    pval_col = "p_combo_t2dobGLP_only", logfc_col = "logFC_combo_t2dobGLP_only",
    s3_subdir = "GLP_only_vs_Neither", file_suffix = "glponly_neither"
  ),
  Both_vs_Neither = list(
    subset_cond = "combo_t2dob %in% c('Both', 'Neither') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "Neither", needs_combo = TRUE,
    pval_col = "p_combo_t2dobBoth", logfc_col = "logFC_combo_t2dobBoth",
    s3_subdir = "Both_vs_Neither", file_suffix = "both_neither"
  ),
  
  # Combo comparisons - Both vs single treatments
  Both_vs_SGLT2i_only = list(
    subset_cond = "combo_t2dob %in% c('Both', 'SGLT2i_only') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "SGLT2i_only", needs_combo = TRUE,
    pval_col = "p_combo_t2dobBoth", logfc_col = "logFC_combo_t2dobBoth",
    s3_subdir = "Both_vs_SGLT2i_only", file_suffix = "both_sglt2ionly"
  ),
  Both_vs_GLP_only = list(
    subset_cond = "combo_t2dob %in% c('Both', 'GLP_only') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "GLP_only", needs_combo = TRUE,
    pval_col = "p_combo_t2dobBoth", logfc_col = "logFC_combo_t2dobBoth",
    s3_subdir = "Both_vs_GLP_only", file_suffix = "both_glponly"
  ),
  
  # Single treatments comparison
  SGLT2i_only_vs_GLP_only = list(
    subset_cond = "combo_t2dob %in% c('SGLT2i_only', 'GLP_only') & !is.na(combo_t2dob)",
    group_var = "combo_t2dob", ref_level = "GLP_only", needs_combo = TRUE,
    pval_col = "p_combo_t2dobSGLT2i_only", logfc_col = "logFC_combo_t2dobSGLT2i_only",
    s3_subdir = "SGLT2i_only_vs_GLP_only", file_suffix = "sglt2ionly_glponly"
  ),
  
  # Combo comparisons within nonDKD (ACR >= 100) - vs HC
  nonDKD_100_Neither_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & combo_t2dob %in% c('Neither', 'HC') & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_Neither_vs_HC", file_suffix = "nondkd100_neither_hc"
  ),
  nonDKD_100_SGLT2i_only_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & combo_t2dob %in% c('SGLT2i_only', 'HC') & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_SGLT2i_only_vs_HC", file_suffix = "nondkd100_sglt2ionly_hc"
  ),
  nonDKD_100_GLP_only_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & combo_t2dob %in% c('GLP_only', 'HC') & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_GLP_only_vs_HC", file_suffix = "nondkd100_glponly_hc"
  ),
  nonDKD_100_Both_vs_HC = list(
    subset_cond = "dkd_group_100_hc != 'DKD' & combo_t2dob %in% c('Both', 'HC') & !is.na(dkd_group_100_hc)",
    group_var = "dkd_group_100_hc", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_dkd_group_100_hcnon_DKD", logfc_col = "logFC_dkd_group_100_hcnon_DKD",
    s3_subdir = "nonDKD_100_Both_vs_HC", file_suffix = "nondkd100_both_hc"
  ),
  
  # Combo comparisons within nonDKD (ACR >= 30) - vs HC
  nonDKD_30_Neither_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & combo_t2dob %in% c('Neither', 'HC') & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_Neither_vs_HC", file_suffix = "nondkd30_neither_hc"
  ),
  nonDKD_30_SGLT2i_only_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & combo_t2dob %in% c('SGLT2i_only', 'HC') & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_SGLT2i_only_vs_HC", file_suffix = "nondkd30_sglt2ionly_hc"
  ),
  nonDKD_30_GLP_only_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & combo_t2dob %in% c('GLP_only', 'HC') & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_GLP_only_vs_HC", file_suffix = "nondkd30_glponly_hc"
  ),
  nonDKD_30_Both_vs_HC = list(
    subset_cond = "dkd_group_30_hc != 'DKD' & combo_t2dob %in% c('Both', 'HC') & !is.na(dkd_group_30_hc)",
    group_var = "dkd_group_30_hc", ref_level = "HC", needs_combo = TRUE,
    pval_col = "p_dkd_group_30_hcnon_DKD", logfc_col = "logFC_dkd_group_30_hcnon_DKD",
    s3_subdir = "nonDKD_30_Both_vs_HC", file_suffix = "nondkd30_both_hc"
  )
)

celltype_groups <- list(
  aPT = "aPT",
  `PT-S1/S2` = "PT-S1/S2",
  `PT-S3` = "PT-S3",
  `PT-1` = "PT-1",
  `PT-2` = "PT-2",
  `PT-3` = "PT-3",
  `PT-4` = "PT-4",
  `PT-5` = "PT-5",
  `C-TAL-1` = "C-TAL-1",
  `C-TAL-2` = "C-TAL-2",
  aTAL = "aTAL",
  dTAL = "dTAL",
  `EC-AVR` = "EC-AVR",
  `EC-GC`  = "EC-GC",
  `EC-PTC` = "EC-PTC",
  `EC-AEA` = "EC-AEA",
  `EC-LYM` = "EC-LYM",
  `EC/VSMC` = "EC/VSMC"
)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


for (folder in names(folders)) {
  config <- analysis_config[[folder]]
  
  for (cell in names(celltype_groups)) {
    
    # Construct S3 input and output keys
    s3_key_raw <- sprintf(
      "Projects/CKD/RH_RH2/Results/nebula/%s/%s/%s_rh_rh2_imp_nebula_kpmp_%s.rds",
      folder, cell, cell, folders[folder]
    )
    
    s3_key_processed <- sprintf(
      "Projects/CKD/RH_RH2/Results/nebula/%s/%s/%s_rh_rh2_imp_nebula_kpmp_%s_processed.rds",
      config$s3_subdir, cell, cell, config$file_suffix
    )
    
    message("\n=====================================")
    message(sprintf("Starting: folder=%s | cell=%s", folder, cell))
    
    # ---- STEP 1: Read raw S3 file ----
    df <- tryCatch({
      s3readRDS(object = s3_key_raw, bucket = "scrna", region = "")
    }, error = function(e) {
      message("❌ FAILED TO READ RAW FILE")
      message(sprintf("Path: %s", s3_key_raw))
      message(sprintf("Details: %s", e$message))
      return(NULL)
    })
    
    if (is.null(df)) {
      message(sprintf("➡️ Skipping: %s | %s\n", folder, cell))
      next
    }
    
    # ---- STEP 2: Process results ----
    processed <- tryCatch({
      process_nebula_results(df,
                             pval_col = config$pval_col,
                             logfc_col = config$logfc_col)
    }, error = function(e) {
      message("❌ FAILED TO PROCESS RESULTS")
      message(sprintf("Details: %s", e$message))
      return(NULL)
    })
    
    if (is.null(processed)) next
    
    # ---- STEP 3: Annotate with gene info ----
    gene_info <- tryCatch({
      getBM(attributes = c("hgnc_symbol", "description", "gene_biotype"),
            filters = "hgnc_symbol",
            values = processed$results$Gene,
            mart = mart) %>%
        dplyr::rename(Gene = hgnc_symbol)
    }, error = function(e) {
      message("⚠️ WARNING: Gene annotation failed, saving unannotated results")
      processed$results  # fallback
    })
    
    annotated_df <- processed$results %>%
      left_join(gene_info, by = "Gene")
    
    # ---- STEP 4: Save processed results to S3 ----
    message(sprintf("Saving processed results → %s", s3_key_processed))
    
    save_ok <- tryCatch({
      s3saveRDS(annotated_df,
                object = s3_key_processed,
                bucket = "scrna",
                region = "",
                multipart = TRUE)
      TRUE
    }, error = function(e) {
      message("❌ FAILED TO SAVE PROCESSED FILE")
      message(sprintf("Path: %s", s3_key_processed))
      message(sprintf("Details: %s", e$message))
      return(FALSE)
    })
    
    if (!save_ok) {
      next
    } else {
      message(sprintf("✅ SUCCESS: %s | %s", folder, cell))
    }
  }
}


