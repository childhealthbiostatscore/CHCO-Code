##### Pierre Analysis 

library(tidyverse)
library(stringr)
library(ggplot2)
library(Seurat)
library(SeuratObject)




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


dat <- dat %>% 
  filter(group %in% c('Type 2 Diabetes', 'Lean Control'))


rh <- data.table::fread("/Users/netio/Documents/UofW/Projects/Pierre_Work/circulating_cell_data/RENALHEIR-PavelSerumCellTypes_DATA_2025-10-02_1324.csv")
improve <- data.table::fread('/Users/netio/Documents/UofW/Projects/Pierre_Work/circulating_cell_data/IMPROVET2D-CirculatingCellData_DATA_2025-10-02_1328.csv')
croc <- data.table::fread('/Users/netio/Documents/UofW/Projects/Pierre_Work/circulating_cell_data/CROCODILE-CirculatingCellData_DATA_2025-10-02_1338.csv')

rh <- rh %>% 
  filter(!is.na(pavel_neutrophil)) %>% 
  dplyr::select(record_id = subject_id, mrn = mr_number, pavel_neutrophil, pavel_lymphocyte, pavel_monocyte)

improve <- improve %>% 
  dplyr::summarize(mrn = mean(mr_number, na.rm=T), 
            pavel_neutrophil = mean(pavel_neutrophil, na.rm=T), 
            pavel_lymphocyte = mean(pavel_lymphocyte, na.rm=T), 
            pavel_monocyte = mean(pavel_monocyte, na.rm=T),
            .by = subject_id) %>% 
  filter(!is.na(pavel_neutrophil)) %>% 
  dplyr::select(record_id = subject_id, mrn, pavel_neutrophil, pavel_lymphocyte, pavel_monocyte)

croc <- croc %>% 
  filter(!is.na(pavel_neutrophils)) %>% 
  dplyr::select(record_id, mrn, pavel_neutrophil = pavel_neutrophils, pavel_lymphocyte, pavel_monocyte)

croc$record_id <- as.character(croc$record_id)

full <- bind_rows(list(rh, improve, croc))



dat_small <- dat %>% filter(mrn %in% full$mrn)


dat_small <- dat_small %>% 
  dplyr::select(mrn, age, hba1c, bmi, total_cholesterol = cholesterol, ldl, hdl, triglycerides)

final <- dat_small %>% left_join(full, by = 'mrn')

names(final) <- c('id', 'age', 'hba1c', 'bmi', 'total_cholesterol', 'ldl', 'hdl', 'triglycerides', 'record', 
                  'neutrophils', 'lymphocytes', 'monocytes')

final <- final %>% dplyr::select(-record)

write.table(final, '/Users/netio/Documents/UofW/Projects/Pierre_Work/circulating_cell_data/UW_participants_formatted.csv', row.names=F, 
            quote=F, sep=',')








############# Cell Type Distributions


group_labels <- data.table::fread("/Users/netio/Documents/UofW/Projects/Pierre_Work/circulating_cell_data/immunotype_results_2025-10-02.csv")
names(group_labels) <- c('mrn', 'cluster')

group_round2 <- readxl::read_excel('/Users/netio/Downloads/CBC_BiopsyParticipants_Final_RHandImprove.xlsx')
group_round2 <- group_round2 %>% filter(!is.na(`Monocytes (ABS #)`)) %>% 
  dplyr::select(record_id = `Subject ID`, mrn = MRN)

group_round2$cluster <- c('MIND', 'SIND', 'MIND', 'LYDD', 'MIND') 

group_round2 <- group_round2 %>% 
  dplyr::select(mrn, cluster)


group_combined <- bind_rows(group_labels, group_round2) %>% 
  filter(!duplicated(mrn))





file_path <- "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/scRNA/data_raw/PB_90_RPCAFix_Metadata_LR_RM_kpmpV1labelled.rds"

so_kpmp_sc <- readRDS(file_path)


dat_small <- dat %>% filter(mrn %in% group_combined$mrn) %>% 
  dplyr::select(mrn, record_id)


so_kpmp_sc <- subset(so_kpmp_sc, subset = record_id %in% dat_small$record_id)


group_combined <- group_combined %>% left_join(dat_small, by='mrn') %>%
  filter(!duplicated(mrn))



#group_combined$record_id[which(group_combined$record_id == 'SWHT_33')] <- 'IT_19'
#group_combined$record_id[which(group_combined$record_id == 'IT_08')] <- 'RH-60-T'
#group_combined$record_id[which(group_combined$record_id == 'IT_10')] <- 'RH-66-T'
#group_combined$record_id[which(group_combined$record_id == 'SWHT_17')] <- 'RH-93-T'


group_mapping <- group_combined


# Calculate total cells per person
total_cells_per_person <- meta.data %>%
  group_by(record_id) %>%
  summarise(total_cells = n(), .groups = 'drop')

# Calculate monocyte counts from KPMP_celltype
monocyte_counts <- meta.data %>% 
  filter(KPMP_celltype == "MON") %>%  # Filter to only monocytes
  group_by(record_id, cluster_mono) %>% 
  summarize(n_monocytes = n(), .groups = 'drop') %>%
  left_join(total_cells_per_person, by = "record_id") %>%
  mutate(pct_monocytes = (n_monocytes / total_cells) * 100)

# View the updated data
head(monocyte_counts)

# Plot monocyte counts
ggplot(monocyte_counts, aes(x = cluster_mono, y = n_monocytes)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Monocyte Counts per Person by Group",
       x = "Group", 
       y = "Number of Monocytes (KPMP_celltype)") +
  theme_minimal()

# Plot monocyte percentages
ggplot(monocyte_counts, aes(x = cluster_mono, y = pct_monocytes)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Monocyte Percentage per Person by Group",
       x = "Group", 
       y = "% Monocytes (of all cells)") +
  theme_minimal()

# Summary statistics
monocyte_summary <- monocyte_counts %>% 
  group_by(cluster_mono) %>% 
  summarize(
    mean_monocytes = mean(n_monocytes, na.rm = T), 
    median_monocytes = median(n_monocytes, na.rm = T), 
    sd_monocytes = sd(n_monocytes, na.rm = T),
    mean_pct_monocytes = mean(pct_monocytes, na.rm = T),
    median_pct_monocytes = median(pct_monocytes, na.rm = T),
    sd_pct_monocytes = sd(pct_monocytes, na.rm = T),
    n_people = n()
  )

print(monocyte_summary)





### Analyses for Pierre 

# =============================================================================
# Analysis of SIND/LYDD Blood Signatures in Kidney Tissue Immune Cells
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(patchwork)

# =============================================================================
# 1. SETUP AND DATA PREPARATION
# =============================================================================

# Load your Seurat object (adjust path as needed)
# seurat_obj <- readRDS("path/to/your/seurat_object.rds")

# Define the 141 DEGs from SIND vs MIND blood analysis
# (You'll need to replace this with your actual list from the Excel file)

bulk_141_genes <- readxl::read_xlsx('/Users/netio/Documents/UofW/Projects/Pierre_Work/intercept_rna_NFS4_SIND_vs_MIND_protein_coding.xlsx',
                                    sheet = 1)
bulk_141_genes <- head(bulk_141_genes$gene_name, 142)

sind_141_genes <- bulk_141_genes

# Define immune cell populations from KPMP atlas
# CORRECT - what you actually have:

lymphoid_cells <- c("pDC", "cDC", "CD4+ T", "CD8+ T", "cycT", "NK", "B")
myeloid_cells <- c("MON", "MAC", "MC")
immune_cells_all <- c(lymphoid_cells, myeloid_cells) 

# =============================================================================
# 2. FILTER TO IMMUNE CELLS ONLY
# =============================================================================

# Filter Seurat object to immune cells
immune_seurat <- subset(so_kpmp_sc, 
                        subset = KPMP_celltype %in% immune_cells_all)

# Add patient group information
immune_seurat$patient_group <- group_mapping$cluster[
  match(immune_seurat$record_id, group_mapping$record_id)
]

# Create immune cell category (lymphoid vs myeloid)
immune_seurat$immune_category <- ifelse(
  immune_seurat$KPMP_celltype %in% lymphoid_cells,
  "Lymphoid", "Myeloid"
)

cat("Immune cells in dataset:\n")
table(immune_seurat$KPMP_celltype)
cat("\nPatient groups:\n")
table(immune_seurat$patient_group)

# =============================================================================
# 3. DOT PLOT: 141 GENES ACROSS IMMUNE CELL TYPES
# =============================================================================

# Check which genes are present in the data
genes_present <- sind_141_genes[sind_141_genes %in% rownames(immune_seurat)]
genes_missing <- sind_141_genes[!sind_141_genes %in% rownames(immune_seurat)]

cat("\nGenes present:", length(genes_present), "/", length(sind_141_genes), "\n")
if(length(genes_missing) > 0) {
  cat("Missing genes:", head(genes_missing, 20), "...\n")
}

# Reorder cell types to group lymphoid and myeloid together
immune_seurat$KPMP_celltype <- factor(immune_seurat$KPMP_celltype,
                                      levels = c(lymphoid_cells, myeloid_cells))

# Create dot plot with cell types ordered by category
p1 <- DotPlot(immune_seurat, 
              features = genes_present,
              group.by = "KPMP_celltype",
              dot.scale = 8,
              cols = c("lightgrey", "red")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 10)) +
  labs(title = "SIND Blood Signature (141 DEGs) in Kidney Immune Cells",
       subtitle = "Lymphoid cells (pDC-B) | Myeloid cells (MON-MC)")

ggsave("dotplot_141genes_immune_cells.pdf", p1, width = 20, height = 8)

# =============================================================================
# 4. UMAP: HIGHLIGHT IMMUNE CELLS
# =============================================================================

# UMAP of all cells with immune cells highlighted
so_kpmp_sc$is_immune <- so_kpmp_sc$KPMP_celltype %in% immune_cells_all

p2 <- DimPlot(so_kpmp_sc, 
              group.by = "is_immune",
              cols = c("grey90", "red"),
              pt.size = 0.5) +
  labs(title = "Immune Cells in Kidney Tissue")

# UMAP of immune cells only, colored by cell type
p3 <- DimPlot(immune_seurat, 
              group.by = "KPMP_celltype",
              label = TRUE,
              repel = TRUE,
              pt.size = 1) +
  labs(title = "Kidney Immune Cell Types")

# UMAP colored by patient group
p4 <- DimPlot(immune_seurat,
              group.by = "patient_group",
              pt.size = 1) +
  labs(title = "Immune Cells by Patient Group (MIND/SIND/LYDD)")

umap_combined <- (p2 | p3) / p4
ggsave("umap_immune_cells_overview.pdf", umap_combined, width = 16, height = 12)

# =============================================================================
# 5. FEATURE PLOTS: TOP EXPRESSED GENES
# =============================================================================

# Find top expressed genes from the 141 list in immune cells
avg_exp <- AverageExpression(immune_seurat, 
                             features = genes_present,
                             group.by = "KPMP_celltype")

# Get top 12 genes by average expression
top_genes <- names(sort(rowMeans(avg_exp$RNA), decreasing = TRUE))[1:12]

p5 <- FeaturePlot(immune_seurat,
                  features = top_genes,
                  ncol = 4,
                  pt.size = 0.5,
                  order = TRUE) &
  theme(legend.position = "right")

ggsave("featureplot_top12_genes_immune.pdf", p5, width = 16, height = 12)

# =============================================================================
# 6. HEATMAP: GENE EXPRESSION BY CELL TYPE AND PATIENT GROUP
# =============================================================================

# Calculate average expression by cell type
avg_by_celltype <- AverageExpression(immune_seurat,
                                     features = genes_present,
                                     group.by = "KPMP_celltype",
                                     slot = "data")

mat1 <- as.matrix(avg_by_celltype$RNA)

# Remove genes with zero variance or all zeros before scaling
gene_vars <- apply(mat1, 1, var)
genes_to_keep <- !is.na(gene_vars) & gene_vars > 0
mat1_filtered <- mat1[genes_to_keep, ]

# Scale the filtered matrix
mat1_scaled <- t(scale(t(mat1_filtered)))

# Remove any remaining NA/Inf values
mat1_scaled[is.na(mat1_scaled)] <- 0
mat1_scaled[is.infinite(mat1_scaled)] <- 0

# Create annotation for lymphoid vs myeloid
cell_anno <- data.frame(
  CellType = colnames(mat1_scaled),
  Category = ifelse(colnames(mat1_scaled) %in% lymphoid_cells, "Lymphoid", "Myeloid")
)
rownames(cell_anno) <- colnames(mat1_scaled)

ha_col <- HeatmapAnnotation(
  Category = cell_anno$Category,
  col = list(Category = c("Lymphoid" = "#4DAF4A", "Myeloid" = "#E41A1C"))
)

pdf("heatmap_141genes_by_celltype.pdf", width = 10, height = 16)
Heatmap(mat1_scaled,
        name = "Scaled\nExpression",
        top_annotation = ha_col,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 10),
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        column_title = paste0("SIND Signature Genes Across Kidney Immune Cells\n(", 
                              nrow(mat1_scaled), " genes shown)"))
dev.off()

cat("\nHeatmap created with", nrow(mat1_scaled), "genes (removed", 
    sum(!genes_to_keep), "genes with zero variance)\n")

# =============================================================================
# 7. STATISTICAL TEST: ENRICHMENT PATTERN
# =============================================================================


test_enrichment <- function(seurat_obj, genes, cell_type_col) {
  results_list <- list()
  cell_types <- unique(seurat_obj[[cell_type_col]][,1])
  
  # Get expression data - compatible with Seurat v5
  expr_data <- GetAssayData(seurat_obj, slot = "data", assay = "RNA")
  
  for(ct in cell_types) {
    # Create cell type label
    seurat_obj$is_celltype <- seurat_obj[[cell_type_col]][,1] == ct
    
    # Test each gene
    ct_results <- data.frame()
    for(gene in genes) {
      if(gene %in% rownames(seurat_obj)) {
        # Extract expression values
        expr_in <- expr_data[gene, seurat_obj$is_celltype]
        expr_out <- expr_data[gene, !seurat_obj$is_celltype]
        
        test <- wilcox.test(expr_in, expr_out)
        
        mean_in <- mean(expr_in)
        mean_out <- mean(expr_out)
        
        ct_results <- rbind(ct_results, data.frame(
          gene = gene,
          cell_type = ct,
          mean_in_celltype = mean_in,
          mean_other_cells = mean_out,
          log2FC = log2((mean_in + 0.01) / (mean_out + 0.01)),
          p_value = test$p.value
        ))
      }
    }
    
    ct_results$padj <- p.adjust(ct_results$p_value, method = "BH")
    results_list[[ct]] <- ct_results
  }
  
  return(bind_rows(results_list))
}

enrichment_results <- test_enrichment(immune_seurat, genes_present, "KPMP_celltype")

# Summary: how many genes are significantly enriched per cell type
enrichment_summary <- enrichment_results %>%
  filter(padj < 0.05, log2FC > 0) %>%
  group_by(cell_type) %>%
  summarise(
    n_enriched_genes = n(),
    prop_enriched = n() / length(genes_present)
  ) %>%
  arrange(desc(n_enriched_genes))

print(enrichment_summary)
write.csv(enrichment_results, "enrichment_test_141genes_immune_cells.csv", row.names = FALSE)

# =============================================================================
# 8. COMPARISON: SIND vs MIND PATIENTS IN TISSUE
# =============================================================================

# Subset to SIND and MIND patients only
sind_mind_immune <- subset(immune_seurat, 
                           subset = patient_group %in% c("SIND", "MIND"))

cat("\nSIND vs MIND comparison:\n")
cat("Total cells:", ncol(sind_mind_immune), "\n")
cat("SIND cells:", sum(sind_mind_immune$patient_group == "SIND"), "\n")
cat("MIND cells:", sum(sind_mind_immune$patient_group == "MIND"), "\n\n")

# Check cell counts per cell type
celltype_counts <- table(sind_mind_immune$KPMP_celltype, sind_mind_immune$patient_group)
print("Cell counts by type and group:")
print(celltype_counts)

# ===================================
# 8A. GROUP-WIDE COMPARISON (ALL IMMUNE CELLS)
# ===================================

cat("\n=== OVERALL IMMUNE CELL COMPARISON ===\n")
Idents(sind_mind_immune) <- "patient_group"

overall_degs <- FindMarkers(sind_mind_immune,
                            ident.1 = "SIND",
                            ident.2 = "MIND",
                            test.use = "wilcox",
                            logfc.threshold = 0.1,
                            min.pct = 0.1)

overall_degs$gene <- rownames(overall_degs)
overall_degs$cell_type <- "All_Immune"
cat("Found", nrow(overall_degs), "DEGs across all immune cells\n")

# Save overall results
write.csv(overall_degs, "tissue_DEGs_SIND_vs_MIND_ALL_immune_cells.csv", row.names = FALSE)

# Check overlap with blood signature
overall_sig_genes <- overall_degs %>%
  filter(p_val_adj < 0.05) %>%
  pull(gene)

overlap_overall <- intersect(genes_present, overall_sig_genes)
cat("\nOverlap with blood SIND signature (ALL immune cells):\n")
cat(length(overlap_overall), "genes out of", length(genes_present), "blood signature genes\n")
if(length(overlap_overall) > 0) {
  print(overlap_overall)
}

# ===================================
# 8B. LYMPHOID vs MYELOID COMPARISONS
# ===================================

cat("\n=== LYMPHOID CELLS COMPARISON ===\n")
sind_mind_lymphoid <- subset(sind_mind_immune, subset = immune_category == "Lymphoid")
if(ncol(sind_mind_lymphoid) > 10) {
  Idents(sind_mind_lymphoid) <- "patient_group"
  
  lymphoid_degs <- FindMarkers(sind_mind_lymphoid,
                               ident.1 = "SIND",
                               ident.2 = "MIND",
                               test.use = "wilcox",
                               logfc.threshold = 0.1,
                               min.pct = 0.1)
  lymphoid_degs$gene <- rownames(lymphoid_degs)
  lymphoid_degs$cell_type <- "Lymphoid"
  cat("Found", nrow(lymphoid_degs), "DEGs in lymphoid cells\n")
  write.csv(lymphoid_degs, "tissue_DEGs_SIND_vs_MIND_LYMPHOID.csv", row.names = FALSE)
  
  overlap_lymphoid <- intersect(genes_present, lymphoid_degs$gene[lymphoid_degs$p_val_adj < 0.05])
  cat("Overlap with blood signature:", length(overlap_lymphoid), "genes\n")
}

cat("\n=== MYELOID CELLS COMPARISON ===\n")
sind_mind_myeloid <- subset(sind_mind_immune, subset = immune_category == "Myeloid")
if(ncol(sind_mind_myeloid) > 10) {
  Idents(sind_mind_myeloid) <- "patient_group"
  
  myeloid_degs <- FindMarkers(sind_mind_myeloid,
                              ident.1 = "SIND",
                              ident.2 = "MIND",
                              test.use = "wilcox",
                              logfc.threshold = 0.1,
                              min.pct = 0.1)
  myeloid_degs$gene <- rownames(myeloid_degs)
  myeloid_degs$cell_type <- "Myeloid"
  cat("Found", nrow(myeloid_degs), "DEGs in myeloid cells\n")
  write.csv(myeloid_degs, "tissue_DEGs_SIND_vs_MIND_MYELOID.csv", row.names = FALSE)
  
  overlap_myeloid <- intersect(genes_present, myeloid_degs$gene[myeloid_degs$p_val_adj < 0.05])
  cat("Overlap with blood signature:", length(overlap_myeloid), "genes\n")
}

# ===================================
# 8C. INDIVIDUAL CELL TYPE COMPARISONS
# ===================================

cat("\n=== INDIVIDUAL CELL TYPE COMPARISONS ===\n")
Idents(sind_mind_immune) <- "patient_group"

tissue_deg_list <- list()
for(ct in immune_cells_all) {
  # Check if cell type exists in filtered data
  if(!(ct %in% sind_mind_immune$KPMP_celltype)) {
    cat("Skipping", ct, "- not present in SIND/MIND samples\n")
    next
  }
  
  ct_cells <- subset(sind_mind_immune, subset = KPMP_celltype == ct)
  
  # Check if we have enough cells in BOTH groups
  n_sind <- sum(ct_cells$patient_group == "SIND")
  n_mind <- sum(ct_cells$patient_group == "MIND")
  
  if(ncol(ct_cells) > 10 && n_sind >= 3 && n_mind >= 3) {
    cat("Testing", ct, "- SIND:", n_sind, "cells, MIND:", n_mind, "cells\n")
    tryCatch({
      degs <- FindMarkers(ct_cells,
                          ident.1 = "SIND",
                          ident.2 = "MIND",
                          test.use = "wilcox",
                          logfc.threshold = 0.1,
                          min.pct = 0.1)
      degs$gene <- rownames(degs)
      degs$cell_type <- ct
      tissue_deg_list[[ct]] <- degs
      cat("  Found", nrow(degs), "DEGs\n")
      
      # Check overlap with blood signature
      ct_overlap <- intersect(genes_present, degs$gene[degs$p_val_adj < 0.05])
      cat("  Overlap with blood signature:", length(ct_overlap), "genes\n")
    }, error = function(e) {
      cat("  Error:", e$message, "\n")
    })
  } else {
    cat("Skipping", ct, "- insufficient cells (SIND:", n_sind, ", MIND:", n_mind, ")\n")
  }
}

if(length(tissue_deg_list) > 0) {
  tissue_degs_all <- bind_rows(tissue_deg_list)
  write.csv(tissue_degs_all, "tissue_DEGs_SIND_vs_MIND_by_celltype.csv", row.names = FALSE)
}

# ===================================
# 8D. SUMMARY COMPARISON
# ===================================

cat("\n=== SUMMARY OF OVERLAPS ===\n")
cat("Blood signature genes tested:", length(genes_present), "\n")
cat("Overall immune cells overlap:", length(overlap_overall), "genes\n")
if(exists("overlap_lymphoid")) cat("Lymphoid cells overlap:", length(overlap_lymphoid), "genes\n")
if(exists("overlap_myeloid")) cat("Myeloid cells overlap:", length(overlap_myeloid), "genes\n")

# Create visualization comparing overlaps
if(exists("lymphoid_degs") && exists("myeloid_degs")) {
  comparison_data <- data.frame(
    Category = c("All Immune", "Lymphoid", "Myeloid"),
    Total_DEGs = c(nrow(overall_degs), nrow(lymphoid_degs), nrow(myeloid_degs)),
    Sig_DEGs = c(sum(overall_degs$p_val_adj < 0.05),
                 sum(lymphoid_degs$p_val_adj < 0.05),
                 sum(myeloid_degs$p_val_adj < 0.05)),
    Blood_Overlap = c(length(overlap_overall),
                      length(overlap_lymphoid),
                      length(overlap_myeloid))
  )
  
  p_comparison <- ggplot(comparison_data, aes(x = Category, y = Blood_Overlap, fill = Category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Blood_Overlap), vjust = -0.5) +
    labs(title = "Blood SIND Signature Genes Found in Kidney Tissue",
         subtitle = "SIND vs MIND comparison",
         x = "Cell Category",
         y = "Number of Overlapping Genes") +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave("overlap_blood_tissue_comparison.pdf", p_comparison, width = 8, height = 6)
}


# Get nominally significant overlaps (p < 0.05, not adjusted)
overlap_nominal_all <- overall_degs %>%
  filter(gene %in% genes_present, p_val < 0.05) %>%
  arrange(p_val)

cat("\n=== NOMINALLY SIGNIFICANT OVERLAPPING GENES ===\n")
cat("Found", nrow(overlap_nominal_all), "genes with nominal p < 0.05\n\n")

# Load full blood data
bulk_141_genes_full <- readxl::read_xlsx('C:/Users/netio/Documents/UofW/Projects/Pierre_Work/intercept_rna_NFS4_SIND_vs_MIND_protein_coding.xlsx',
                                         sheet = 1)

if(nrow(overlap_nominal_all) > 0) {
  print(overlap_nominal_all %>% 
          select(gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2))
  
  # Get blood data for these genes to compare directions
  blood_data <- bulk_141_genes_full %>%
    filter(gene_name %in% overlap_nominal_all$gene)
  
  # Check direction concordance - FIXED column name
  comparison_table <- overlap_nominal_all %>%
    select(gene, tissue_log2FC = avg_log2FC, tissue_pval = p_val, tissue_padj = p_val_adj) %>%
    left_join(blood_data %>% select(gene = gene_name, blood_log2FC = log2FoldChange, blood_padj = padj),
              by = "gene") %>%
    mutate(same_direction = sign(tissue_log2FC) == sign(blood_log2FC))
  
  print("\nDirection Concordance:")
  print(comparison_table)
  
  # Visualize these genes
  if(nrow(overlap_nominal_all) <= 10) {
    # Violin plots for each gene
    p_genes <- VlnPlot(sind_mind_immune,
                       features = overlap_nominal_all$gene,
                       group.by = "patient_group",
                       split.by = "patient_group",
                       pt.size = 0.1,
                       ncol = 2) 
    
    ggsave("nominal_overlap_genes_violinplot.pdf", p_genes, width = 10, height = 8)
    
    # Feature plots
    p_features <- FeaturePlot(sind_mind_immune,
                              features = overlap_nominal_all$gene,
                              split.by = "patient_group",
                              ncol = 2)
    
    ggsave("nominal_overlap_genes_featureplot.pdf", p_features, width = 12, height = 8)
    
    # Dot plot across cell types
    p_dot <- DotPlot(sind_mind_immune,
                     features = overlap_nominal_all$gene,
                     group.by = "KPMP_celltype",
                     split.by = "patient_group") +
      RotatedAxis() +
      labs(title = "Nominally Significant Overlapping Genes\nAcross Cell Types")
    
    ggsave("nominal_overlap_genes_by_celltype.pdf", p_dot, width = 10, height = 6)
  }
}

# Statistical note for interpretation
cat("\n=== STATISTICAL NOTE ===\n")
cat("Expected by chance at p < 0.05:", round(length(genes_present) * 0.05, 1), "genes\n")
cat("Observed:", nrow(overlap_nominal_all), "genes\n")
if(nrow(overlap_nominal_all) > length(genes_present) * 0.05) {
  cat("This is MORE than expected by chance - suggests real signal\n")
} else {
  cat("This is consistent with chance - may not represent real signal\n")
}

# Binomial test: is this more than chance?
binom_test <- binom.test(x = nrow(overlap_nominal_all),
                         n = length(genes_present),
                         p = 0.05,
                         alternative = "greater")
cat("\nBinomial test p-value:", binom_test$p.value, "\n")
if(binom_test$p.value < 0.05) {
  cat("Significantly more overlap than expected by chance!\n")
} else {
  cat("Not significantly more than random chance.\n")
}



# =============================================================================
# 9. SUMMARY VISUALIZATION
# =============================================================================

# Create summary barplot showing enrichment across cell types
p6 <- ggplot(enrichment_summary, 
             aes(x = reorder(cell_type, n_enriched_genes), 
                 y = n_enriched_genes,
                 fill = ifelse(cell_type %in% myeloid_cells, "Myeloid", "Lymphoid"))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Number of SIND Blood Signature Genes\nEnriched in Each Kidney Immune Cell Type",
       x = "Cell Type",
       y = "Number of Enriched Genes (padj < 0.05)",
       fill = "Category") +
  scale_fill_manual(values = c("Lymphoid" = "#4DAF4A", "Myeloid" = "#E41A1C")) +
  theme_minimal() +
  theme(text = element_text(size = 12))

ggsave("summary_enrichment_barplot.pdf", p6, width = 8, height = 6)

cat("\n=============================================================================\n")
cat("Analysis complete! Generated files:\n")
cat("1. dotplot_141genes_immune_cells.pdf\n")
cat("2. umap_immune_cells_overview.pdf\n")
cat("3. featureplot_top12_genes_immune.pdf\n")
cat("4. heatmap_141genes_by_celltype.pdf\n")
cat("5. enrichment_test_141genes_immune_cells.csv\n")
cat("6. tissue_DEGs_SIND_vs_MIND_immune_cells.csv\n")
cat("7. summary_enrichment_barplot.pdf\n")
cat("=============================================================================\n")




# =============================================================================
# ADDITIONAL ANALYSES FOR SIND GRANT - STRENGTHENING PRELIMINARY DATA
# =============================================================================
# This script performs additional analyses to support the blood-tissue 
# connection for the SIND signature in kidney immune cells
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(fgsea)
library(pheatmap)
library(ggpubr)
library(tibble)

# Assumes you have already run the main analysis script and have:
# - immune_seurat: Seurat object with immune cells
# - sind_mind_immune: Subset to SIND and MIND patients
# - overall_degs: DEGs from all immune cells SIND vs MIND
# - genes_present: 141 SIND signature genes present in data
# - bulk_141_genes_full: Full blood signature data with log2FC

# =============================================================================
# 1. GENE SET ENRICHMENT ANALYSIS (GSEA)
# =============================================================================

cat("\n=== 1. GENE SET ENRICHMENT ANALYSIS ===\n")

# Prepare ranked gene list from tissue (all immune cells, SIND vs MIND)
gene_ranks <- overall_degs %>%
  mutate(rank_metric = -log10(p_val + 1e-10) * sign(avg_log2FC)) %>%
  arrange(desc(rank_metric)) %>%
  pull(rank_metric, name = gene)

# Your SIND signature as a gene set
sind_signature <- list(SIND_Blood_Signature = genes_present)

# Run GSEA
set.seed(42)
gsea_results <- fgsea(pathways = sind_signature,
                      stats = gene_ranks,
                      minSize = 10,
                      maxSize = 500,
                      nperm = 10000)

cat("\nGSEA Results:\n")
print(gsea_results)

# Plot enrichment
pdf("gsea_sind_signature_kidney.pdf", width = 8, height = 6)
plotEnrichment(sind_signature$SIND_Blood_Signature, gene_ranks) +
  labs(title = "SIND Blood Signature Enrichment in Kidney Tissue",
       subtitle = paste0("NES = ", round(gsea_results$NES, 2), 
                         ", p = ", format.pval(gsea_results$pval, digits = 3)))
dev.off()

cat("GSEA plot saved: gsea_sind_signature_kidney.pdf\n")

# =============================================================================
# 2. BLOOD-TISSUE COMPARISON HEATMAP
# =============================================================================

cat("\n=== 2. BLOOD-TISSUE COMPARISON HEATMAP ===\n")

# Get blood log2FC for the 141 genes
blood_fc <- bulk_141_genes_full %>%
  filter(gene_name %in% genes_present) %>%
  select(gene = gene_name, blood_log2FC = log2FoldChange) %>%
  arrange(desc(blood_log2FC))

# Get tissue log2FC from overall immune comparison
tissue_fc <- overall_degs %>%
  filter(gene %in% genes_present) %>%
  select(gene, tissue_log2FC = avg_log2FC)

# Combine
comparison_matrix <- blood_fc %>%
  left_join(tissue_fc, by = "gene") %>%
  mutate(tissue_log2FC = ifelse(is.na(tissue_log2FC), 0, tissue_log2FC)) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Cap extreme values for better visualization
comparison_matrix[comparison_matrix > 2] <- 2
comparison_matrix[comparison_matrix < -2] <- -2

# Plot
pheatmap(comparison_matrix,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-2, 2, length.out = 101),
         main = "SIND Signature: Blood vs Kidney Tissue Log2FC",
         angle_col = 0,
         fontsize_row = 4,
         filename = "blood_tissue_comparison_heatmap.pdf",
         width = 5,
         height = 12)

cat("Heatmap saved: blood_tissue_comparison_heatmap.pdf\n")

# =============================================================================
# 3. CORRELATION BETWEEN BLOOD AND TISSUE FOLD CHANGES
# =============================================================================

cat("\n=== 3. CORRELATION ANALYSIS ===\n")

correlation_data <- blood_fc %>%
  left_join(tissue_fc, by = "gene") %>%
  filter(!is.na(tissue_log2FC))

# Pearson correlation
cor_test <- cor.test(correlation_data$blood_log2FC, 
                     correlation_data$tissue_log2FC)

cat("\nPearson correlation: r =", round(cor_test$estimate, 3), 
    ", p =", format.pval(cor_test$p.value, digits = 3), "\n")

# Spearman correlation (rank-based, more robust)
cor_test_spearman <- cor.test(correlation_data$blood_log2FC,
                              correlation_data$tissue_log2FC,
                              method = "spearman")

cat("Spearman correlation: rho =", round(cor_test_spearman$estimate, 3),
    ", p =", format.pval(cor_test_spearman$p.value, digits = 3), "\n")

# Plot
p_correlation <- ggplot(correlation_data, 
                        aes(x = blood_log2FC, y = tissue_log2FC)) +
  geom_point(alpha = 0.6, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(title = "Blood vs Tissue Log2 Fold Changes",
       subtitle = paste0("Pearson r = ", round(cor_test$estimate, 3),
                         ", p = ", format.pval(cor_test$p.value, digits = 2)),
       x = "Blood Log2FC (SIND vs MIND)",
       y = "Kidney Tissue Log2FC (SIND vs MIND)") +
  theme_minimal() +
  theme(text = element_text(size = 12))

ggsave("blood_tissue_correlation.pdf", p_correlation, width = 8, height = 6)

# Direction concordance
concordance <- correlation_data %>%
  mutate(same_direction = sign(blood_log2FC) == sign(tissue_log2FC),
         both_upregulated = blood_log2FC > 0 & tissue_log2FC > 0,
         both_downregulated = blood_log2FC < 0 & tissue_log2FC < 0)

cat("\nDirection concordance:\n")
cat("Same direction:", sum(concordance$same_direction), "/", nrow(concordance), 
    "(", round(100*mean(concordance$same_direction), 1), "%)\n")
cat("Both upregulated:", sum(concordance$both_upregulated), "\n")
cat("Both downregulated:", sum(concordance$both_downregulated), "\n")

# Save concordance data
write.csv(concordance, "blood_tissue_concordance.csv", row.names = FALSE)

# =============================================================================
# 4. CELL TYPE SPECIFICITY OF SIND SIGNATURE
# =============================================================================

cat("\n=== 4. CELL TYPE SPECIFICITY ===\n")

# Calculate average expression of 141 genes per cell type
avg_expression_by_celltype <- AverageExpression(immune_seurat,
                                                features = genes_present,
                                                group.by = "KPMP_celltype")

# Get mean expression per cell type
celltype_signature_score <- colMeans(avg_expression_by_celltype$RNA)

# Plot
signature_data <- data.frame(
  CellType = names(celltype_signature_score),
  MeanExpression = celltype_signature_score,
  Category = ifelse(names(celltype_signature_score) %in% myeloid_cells, 
                    "Myeloid", "Lymphoid")
)

p_celltype_expr <- ggplot(signature_data, 
                          aes(x = reorder(CellType, MeanExpression), 
                              y = MeanExpression, 
                              fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "SIND Blood Signature Expression in Kidney Immune Cells",
       subtitle = paste0("Mean expression of ", length(genes_present), " signature genes"),
       x = "Cell Type",
       y = "Mean Expression") +
  scale_fill_manual(values = c("Lymphoid" = "#4DAF4A", "Myeloid" = "#E41A1C")) +
  theme_minimal() +
  theme(text = element_text(size = 12))

ggsave("sind_signature_by_celltype.pdf", p_celltype_expr, width = 8, height = 6)

cat("Highest expression in:", signature_data$CellType[which.max(signature_data$MeanExpression)], "\n")

# =============================================================================
# 5. MODULE SCORE ANALYSIS
# =============================================================================

cat("\n=== 5. MODULE SCORE ANALYSIS ===\n")

# Add module score for SIND signature
immune_seurat <- AddModuleScore(immune_seurat,
                                features = list(genes_present),
                                name = "SIND_Signature")

# Rename the column (AddModuleScore adds a number)
score_col <- grep("SIND_Signature", colnames(immune_seurat@meta.data), value = TRUE)
immune_seurat@meta.data$SIND_Score <- immune_seurat@meta.data[[score_col]]

# Compare scores between SIND and MIND
score_comparison <- immune_seurat@meta.data %>%
  filter(patient_group %in% c("SIND", "MIND")) %>%
  select(patient_group, KPMP_celltype, SIND_Score, immune_category)

# Overall comparison
wilcox_overall <- wilcox.test(SIND_Score ~ patient_group, data = score_comparison)
cat("\nOverall SIND score comparison:\n")
cat("Wilcoxon p-value:", format.pval(wilcox_overall$p.value, digits = 3), "\n")
cat("Mean SIND score:", mean(score_comparison$SIND_Score[score_comparison$patient_group == "SIND"]), "\n")
cat("Mean MIND score:", mean(score_comparison$SIND_Score[score_comparison$patient_group == "MIND"]), "\n")

p_score_overall <- ggplot(score_comparison, 
                          aes(x = patient_group, y = SIND_Score, fill = patient_group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.2, size = 0.5) +
  stat_compare_means(method = "wilcox.test", label.y.npc = 0.95) +
  labs(title = "SIND Blood Signature Score in Kidney Tissue",
       subtitle = "All Immune Cells Combined",
       x = "Patient Group",
       y = "SIND Signature Score") +
  scale_fill_manual(values = c("SIND" = "#E41A1C", "MIND" = "#377EB8")) +
  theme_minimal() +
  theme(text = element_text(size = 12), legend.position = "none")

ggsave("sind_score_overall.pdf", p_score_overall, width = 6, height = 6)

# By cell type
p_score_celltype <- ggplot(score_comparison, 
                           aes(x = KPMP_celltype, y = SIND_Score, fill = patient_group)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_grid(immune_category ~ ., scales = "free_y", space = "free_y") +
  labs(title = "SIND Signature Score by Cell Type",
       x = "Cell Type",
       y = "SIND Signature Score") +
  scale_fill_manual(values = c("SIND" = "#E41A1C", "MIND" = "#377EB8")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("sind_score_by_celltype.pdf", p_score_celltype, width = 10, height = 8)

# Test by cell type
score_by_celltype <- score_comparison %>%
  group_by(KPMP_celltype) %>%
  summarise(
    n_sind = sum(patient_group == "SIND"),
    n_mind = sum(patient_group == "MIND"),
    mean_sind = mean(SIND_Score[patient_group == "SIND"]),
    mean_mind = mean(SIND_Score[patient_group == "MIND"]),
    diff = mean_sind - mean_mind
  ) %>%
  arrange(desc(diff))

cat("\nScore differences by cell type (SIND - MIND):\n")
print(score_by_celltype)

write.csv(score_by_celltype, "sind_score_by_celltype_summary.csv", row.names = FALSE)

# =============================================================================
# 6. TOP 30 BLOOD GENES IN TISSUE
# =============================================================================

cat("\n=== 6. TOP 30 BLOOD GENES ANALYSIS ===\n")

# Get top 30 genes from blood by adjusted p-value
top_blood_genes <- bulk_141_genes_full %>%
  arrange(padj) %>%
  head(30) %>%
  pull(gene_name)

top_blood_genes <- top_blood_genes[top_blood_genes %in% rownames(immune_seurat)]

cat("Testing top", length(top_blood_genes), "genes from blood\n")

# Dot plot for top genes only
p_top_genes <- DotPlot(immune_seurat,
                       features = top_blood_genes,
                       group.by = "KPMP_celltype",
                       dot.scale = 8) +
  RotatedAxis() +
  labs(title = "Top 30 SIND Blood Signature Genes in Kidney Immune Cells",
       subtitle = "Ranked by blood adjusted p-value") +
  theme(axis.text.x = element_text(size = 8))

ggsave("top30_blood_genes_kidney.pdf", p_top_genes, width = 14, height = 6)

# Check overlap for top genes only
overlap_top <- overall_degs %>%
  filter(gene %in% top_blood_genes, p_val < 0.05)

cat("\nTop 30 blood genes with nominal significance (p<0.05) in tissue:", nrow(overlap_top), "\n")
if(nrow(overlap_top) > 0) {
  print(overlap_top %>% select(gene, avg_log2FC, p_val, p_val_adj))
}

# Module score for top 30 only
immune_seurat <- AddModuleScore(immune_seurat,
                                features = list(top_blood_genes),
                                name = "SIND_Top30")

score_col_top30 <- grep("SIND_Top30", colnames(immune_seurat@meta.data), value = TRUE)
immune_seurat@meta.data$SIND_Top30_Score <- immune_seurat@meta.data[[score_col_top30]]

score_comparison_top30 <- immune_seurat@meta.data %>%
  filter(patient_group %in% c("SIND", "MIND"))

wilcox_top30 <- wilcox.test(SIND_Top30_Score ~ patient_group, 
                            data = score_comparison_top30)

cat("\nTop 30 genes signature score:\n")
cat("Wilcoxon p-value:", format.pval(wilcox_top30$p.value, digits = 3), "\n")

# =============================================================================
# 7. EXPRESSION PATTERNS: ARE SIGNATURE GENES EXPRESSED?
# =============================================================================

cat("\n=== 7. GENE EXPRESSION PATTERNS ===\n")

# Check what percentage of signature genes are expressed in tissue
expr_data <- GetAssayData(immune_seurat, slot = "data", assay = "RNA")

genes_expressed <- sapply(genes_present, function(g) {
  if(g %in% rownames(expr_data)) {
    sum(expr_data[g, ] > 0)
  } else {
    0
  }
})

genes_expressed_df <- data.frame(
  gene = names(genes_expressed),
  n_cells_expressed = genes_expressed,
  pct_cells_expressed = 100 * genes_expressed / ncol(immune_seurat)
) %>%
  arrange(desc(n_cells_expressed))

cat("\nExpression summary:\n")
cat("Genes expressed in >10% of cells:", sum(genes_expressed_df$pct_cells_expressed > 10), "/", 
    length(genes_present), "\n")
cat("Genes expressed in >50% of cells:", sum(genes_expressed_df$pct_cells_expressed > 50), "/",
    length(genes_present), "\n")
cat("Genes not expressed:", sum(genes_expressed_df$n_cells_expressed == 0), "/",
    length(genes_present), "\n")

write.csv(genes_expressed_df, "sind_signature_expression_summary.csv", row.names = FALSE)

# Plot expression distribution
p_expression_dist <- ggplot(genes_expressed_df, 
                            aes(x = pct_cells_expressed)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  labs(title = "SIND Signature Gene Expression in Kidney Immune Cells",
       x = "% of Cells Expressing Gene",
       y = "Number of Genes") +
  theme_minimal()

ggsave("signature_gene_expression_distribution.pdf", p_expression_dist, width = 8, height = 6)

# =============================================================================
# 8. SUMMARY FIGURE: MULTI-PANEL FOR GRANT
# =============================================================================

cat("\n=== 8. CREATING SUMMARY FIGURE ===\n")

library(patchwork)

# Panel A: Correlation
panel_a <- p_correlation + 
  labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold", size = 16))

# Panel B: Module score
panel_b <- p_score_overall + 
  labs(tag = "B") +
  theme(plot.tag = element_text(face = "bold", size = 16))

# Panel C: Cell type specificity
panel_c <- p_celltype_expr + 
  labs(tag = "C") +
  theme(plot.tag = element_text(face = "bold", size = 16))

# Panel D: Expression distribution
panel_d <- p_expression_dist +
  labs(tag = "D") +
  theme(plot.tag = element_text(face = "bold", size = 16))

# Combine
summary_figure <- (panel_a | panel_b) / (panel_c | panel_d)

ggsave("summary_figure_for_grant.pdf", summary_figure, width = 14, height = 12)
ggsave("summary_figure_for_grant.png", summary_figure, width = 14, height = 12, dpi = 300)

cat("\nSummary figure saved: summary_figure_for_grant.pdf and .png\n")



# =============================================================================
# 8. SUMMARY FIGURE: MULTI-PANEL FOR GRANT (UPDATED)
# =============================================================================

cat("\n=== 8. CREATING SUMMARY FIGURE ===\n")

library(patchwork)

# Panel A: SIND Score by Cell Type (REPLACES correlation plot)
# Improved Panel A with side-by-side facets
p_score_celltype_improved <- ggplot(score_comparison, 
                                    aes(x = KPMP_celltype, y = SIND_Score, fill = patient_group)) +
  geom_boxplot(outlier.size = 0.8, linewidth = 0.6, alpha = 0.8) +
  facet_grid(. ~ immune_category, scales = "free_x", space = "free_x") +  # Changed to side-by-side
  labs(title = "SIND Signature Score by Cell Type",
       x = "Cell Type",
       y = "SIND Signature Score",
       fill = "Patient Group") +
  scale_fill_manual(values = c("SIND" = "#E41A1C", "MIND" = "#377EB8")) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 13, face = "bold"),
        strip.text = element_text(size = 12, face = "bold", color = "white"),
        strip.background = element_rect(fill = c("Lymphoid" = "#4DAF4A", "Myeloid" = "#E41A1C")),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(1, "lines"))

panel_a <- p_score_celltype_improved + 
  labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold", size = 16))

# Panel B: Module score overall
panel_b <- p_score_overall + 
  labs(tag = "B") +
  theme(plot.tag = element_text(face = "bold", size = 16),
        text = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        plot.title = element_text(size = 14, face = "bold"))

# Panel C: Cell type specificity
panel_c <- p_celltype_expr + 
  labs(tag = "C") +
  theme(plot.tag = element_text(face = "bold", size = 16),
        text = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))

# Panel D: Expression distribution
panel_d <- p_expression_dist +
  labs(tag = "D") +
  theme(plot.tag = element_text(face = "bold", size = 16),
        text = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        plot.title = element_text(size = 14, face = "bold"))

# Combine with better layout
summary_figure <- (panel_a | panel_b) / (panel_c | panel_d)

# Save high quality versions
ggsave("summary_figure_for_grant.pdf", 
       summary_figure, 
       width = 16, 
       height = 14,
       dpi = 300,
       device = cairo_pdf)  # Better quality PDF

ggsave("summary_figure_for_grant.png", 
       summary_figure, 
       width = 16, 
       height = 14, 
       dpi = 600,  # High resolution for publication
       bg = "white")

ggsave("summary_figure_for_grant.tiff", 
       summary_figure, 
       width = 16, 
       height = 14, 
       dpi = 600,  # TIFF format often required for journals
       compression = "lzw",
       bg = "white")


















