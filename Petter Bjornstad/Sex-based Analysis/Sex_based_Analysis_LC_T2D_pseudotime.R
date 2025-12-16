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





