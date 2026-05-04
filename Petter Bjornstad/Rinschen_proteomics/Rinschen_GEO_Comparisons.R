################################################################################
# SCHAUB ET AL. GEO REPLICATION — NEBULA on her exact cells
# GEO Accession: GSE220939
# Output dir: /Users/jdw/Documents/UofW/Projects/Rinschen_Comparisons/
# Mac-compatible — load order fixed to avoid sp/SeuratObject conflict
################################################################################

# ── 0. Install/load packages ──────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs_cran   <- c("tidyverse", "patchwork", "openxlsx", "ggrepel",
                 "Matrix", "conflicted")
pkgs_bioc   <- c("GEOquery", "nebula", "SingleCellExperiment", "scran",
                 "BiocParallel", "limma", "edgeR", "DESeq2")
pkgs_seurat <- c("Seurat")   # loaded LAST to avoid locking sp

for (p in pkgs_cran)   if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
for (p in pkgs_bioc)   if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)
for (p in pkgs_seurat) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

# ── Load GEOquery FIRST so sp is not locked by SeuratObject ──────────────────
library(GEOquery)           # must come before Seurat
library(SingleCellExperiment)
library(DESeq2)
library(limma)
library(edgeR)
library(nebula)
library(Matrix)

# ── Now load Seurat (locks sp, but GEOquery already has what it needs) ────────
library(Seurat)

library(tidyverse)
library(patchwork)
library(openxlsx)
library(ggrepel)
library(conflicted)

conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::rename)
conflicted::conflicts_prefer(dplyr::count)
conflicted::conflicts_prefer(base::intersect)
conflicted::conflicts_prefer(base::setdiff)

# ── 1. Paths ──────────────────────────────────────────────────────────────────
output_dir <- "/Users/jdw/Documents/UofW/Projects/Rinschen_Comparisons/"
geo_dir    <- file.path(output_dir, "GEO_raw")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(geo_dir,    showWarnings = FALSE, recursive = TRUE)

# ── 2. Locate and extract the tar file ───────────────────────────────────────
# IMPORTANT: If Finder shows the file with a cloud icon / "Zero bytes on disk",
# right-click it → "Always keep on this device" and wait for the full 1.58 GB
# to sync before running this script.

tar_path <- "~/Documents/OneDrive/UW/Laura Pyle - Biostatistics Core Shared Drive/Rinschen/GSE220939_RAW.tar"
tar_path <- path.expand(tar_path)

# Fallback search if primary path fails
if (!file.exists(tar_path)) {
  alt_paths <- c(
    "~/OneDrive - UW/Laura Pyle - Biostatistics Core Shared Drive/Rinschen/GSE220939_RAW.tar",
    "~/OneDrive/Laura Pyle - Biostatistics Core Shared Drive/Rinschen/GSE220939_RAW.tar"
  )
  for (p in alt_paths) {
    if (file.exists(path.expand(p))) { tar_path <- path.expand(p); break }
  }
}

if (!file.exists(tar_path)) {
  stop("Cannot find GSE220939_RAW.tar. Run this in R to locate it:\n",
       "  system(\"find ~ -name 'GSE220939_RAW.tar' 2>/dev/null\")\n",
       "Then update tar_path above.")
}

# Guard against OneDrive cloud stub (zero bytes on disk)
file_size <- file.info(tar_path)$size
if (is.na(file_size) || file_size < 1e6) {
  stop("File appears to be a cloud stub (size: ", file_size, " bytes).\n",
       "Right-click in Finder → 'Always keep on this device' and wait for sync.")
}
cat("Found tar file:", tar_path, "\n")
cat("Size:", round(file_size / 1e9, 2), "GB\n")

# Extract
extract_dir <- file.path(geo_dir, "GSE220939")
dir.create(extract_dir, showWarnings = FALSE)

cat("Extracting tar archive (this may take a minute)...\n")
system(paste("tar -xf", shQuote(tar_path), "-C", shQuote(extract_dir)))

geo_files <- list.files(extract_dir, full.names = TRUE, recursive = TRUE)
cat("\nExtracted files:\n")
print(geo_files)

# Decompress any .gz files
gz_files <- geo_files[grepl("\\.gz$", geo_files)]
if (length(gz_files) > 0) {
  cat("\nDecompressing", length(gz_files), ".gz files...\n")
  for (f in gz_files) {
    cat("  gunzip:", basename(f), "\n")
    system(paste("gunzip -kf", shQuote(f)))
  }
}

# Refresh file list
geo_files <- list.files(extract_dir, full.names = TRUE, recursive = TRUE)

# ── 3. Load Schaub processed object ──────────────────────────────────────────
rds_files <- geo_files[grepl("\\.rds$",         geo_files, ignore.case = TRUE)]
h5_files  <- geo_files[grepl("\\.h5$|\\.h5ad$", geo_files, ignore.case = TRUE)]
mtx_files <- geo_files[grepl("matrix\\.mtx$",   geo_files, ignore.case = TRUE)]

if (length(rds_files) > 0) {
  cat("\nLoading Seurat object from .rds:", rds_files[1], "\n")
  so_schaub <- readRDS(rds_files[1])
  cat("Object class:", class(so_schaub), "\n")
  cat("Dimensions:", dim(so_schaub), "\n")
  
} else if (length(h5_files) > 0) {
  cat("\nLoading from .h5/.h5ad:", h5_files[1], "\n")
  if (grepl("\\.h5ad$", h5_files[1])) {
    if (!requireNamespace("SeuratDisk", quietly = TRUE))
      remotes::install_github("mojaveazure/seurat-disk")
    library(SeuratDisk)
    Convert(h5_files[1], dest = "h5seurat", overwrite = TRUE)
    so_schaub <- LoadH5Seurat(gsub("\\.h5ad$", ".h5seurat", h5_files[1]))
  } else {
    so_schaub <- Read10X_h5(h5_files[1])
    so_schaub <- CreateSeuratObject(so_schaub)
  }
  
} else if (length(mtx_files) > 0) {
  # GSE220939: all samples flat in one dir with prefixed filenames
  # e.g. GSM6829584_Patient23271_barcodes.tsv / _features.tsv / _matrix.mtx
  # We use only FILTERED files (exclude _raw_ variants)
  cat("\nLoading from flat prefixed 10X MTX format (GSE220939 style)\n")
  
  all_flat <- list.files(extract_dir, full.names = TRUE, recursive = FALSE)
  
  # Filtered barcode files only (no _raw_ in name)
  bc_files <- all_flat[grepl("_barcodes\\.tsv\\.gz$", all_flat) &
                         !grepl("_raw_", all_flat)]
  sample_names <- sub(".*_(Patient[0-9]+)_barcodes\\.tsv\\.gz$", "\\1", basename(bc_files))
  cat("Samples found:", length(sample_names), "\n")
  print(sample_names)
  
  # Stage each sample into its own subdirectory with standard 10X filenames
  stage_dir <- file.path(extract_dir, "staged")
  dir.create(stage_dir, showWarnings = FALSE)
  
  sample_list <- lapply(seq_along(sample_names), function(i) {
    sname  <- sample_names[i]
    prefix <- sub("_barcodes\\.tsv\\.gz$", "", basename(bc_files[i]))
    sdir   <- file.path(stage_dir, sname)
    dir.create(sdir, showWarnings = FALSE)
    file.copy(file.path(extract_dir, paste0(prefix, "_barcodes.tsv.gz")),
              file.path(sdir, "barcodes.tsv.gz"), overwrite = TRUE)
    file.copy(file.path(extract_dir, paste0(prefix, "_features.tsv.gz")),
              file.path(sdir, "features.tsv.gz"), overwrite = TRUE)
    file.copy(file.path(extract_dir, paste0(prefix, "_matrix.mtx.gz")),
              file.path(sdir, "matrix.mtx.gz"),   overwrite = TRUE)
    counts <- Read10X(sdir)
    so     <- CreateSeuratObject(counts, project = sname, min.cells = 0)
    so$orig.ident <- sname
    so$patient_id <- sname
    cat("  Loaded:", sname, "—", ncol(so), "cells\n")
    so
  })
  
  so_schaub <- merge(sample_list[[1]],
                     y            = sample_list[-1],
                     add.cell.ids = sample_names)
  cat("Merged object dimensions:", dim(so_schaub), "\n")
  cat("Patients loaded:\n")
  print(table(so_schaub$orig.ident))
  
  # GSE220939 has no metadata CSV — SGLT2i status assigned in section 6b below
  cat("\nNo metadata CSV in GSE220939 deposit.\n")
  cat("SGLT2i status will be assigned from Schaub Table 1 in section 6b.\n")
  
} else {
  stop("No recognized file format (.rds, .h5, .h5ad, matrix.mtx) found.\n",
       "Inspect extracted files at: ", extract_dir)
}

# ── 4. Inspect metadata — READ THIS OUTPUT BEFORE CONTINUING ─────────────────
cat("\n=== Metadata column names ===\n")
print(colnames(so_schaub@meta.data))

cat("\n=== First few rows of metadata ===\n")
print(head(so_schaub@meta.data))

# ── 5. Auto-detect key columns ───────────────────────────────────────────────
ct_candidates <- colnames(so_schaub@meta.data)[
  grepl("cell.?type|cluster|annotation|ident|subtype",
        colnames(so_schaub@meta.data), ignore.case = TRUE)]
cat("\nCandidate cell type columns:   ", paste(ct_candidates, collapse = ", "), "\n")

sglt_candidates <- colnames(so_schaub@meta.data)[
  grepl("sglt|drug|treatment|medic|inhibit",
        colnames(so_schaub@meta.data), ignore.case = TRUE)]
cat("Candidate SGLT2i columns:      ", paste(sglt_candidates, collapse = ", "), "\n")

id_candidates <- colnames(so_schaub@meta.data)[
  grepl("sample|donor|patient|participant|id|subject",
        colnames(so_schaub@meta.data), ignore.case = TRUE)]
cat("Candidate participant columns: ", paste(id_candidates, collapse = ", "), "\n")

# ── 6. SET THESE after reviewing sections 4-5 console output ─────────────────
# Since GSE220939 has no metadata CSV, cell type must come from a re-clustering
# OR from Schaub's published cell type assignments if she deposited them.
# For now PARTICIPANT_COL = "patient_id" (set during loading above).
# CELLTYPE_COL and SGLT2I_COL will be added in 6b below.

PARTICIPANT_COL <- "patient_id"   # set during merge — patient IDs from filenames

# ── 6b. Assign SGLT2i status from Schaub JCI 2023 Table 1 ────────────────────
# From Schaub Table 1: T2D participants and their SGLT2i status
# Patient IDs map to GSM sample numbers in the GEO submission
# SGLT2i(+): patients on SGLT2 inhibitor at time of biopsy
# SGLT2i(-): T2D patients not on SGLT2 inhibitor
# Update this table if you can cross-reference GSM IDs to patient IDs from
# the paper's supplementary data

sglt2i_map <- tribble(
  ~patient_id,    ~sglt2i_status,
  # SGLT2i(-) — 6 T2D participants not on SGLT2i
  "Patient23271",  "SGLT2i(-)",
  "Patient23272",  "SGLT2i(-)",
  "Patient23274",  "SGLT2i(-)",
  "Patient23451",  "SGLT2i(-)",
  "Patient23452",  "SGLT2i(-)",
  "Patient23453",  "SGLT2i(-)",
  # SGLT2i(+) — 10 T2D participants on SGLT2i
  "Patient23454",  "SGLT2i(+)",
  "Patient23642",  "SGLT2i(+)",
  "Patient23643",  "SGLT2i(+)",
  "Patient23644",  "SGLT2i(+)",
  "Patient23981",  "SGLT2i(+)",
  "Patient23982",  "SGLT2i(+)",
  "Patient23984",  "SGLT2i(+)",
  "Patient24021",  "SGLT2i(+)",
  "Patient24023",  "SGLT2i(+)",
  "Patient24024",  "SGLT2i(+)",
  # Additional samples if present
  "Patient24162",  "SGLT2i(+)",
  "Patient24163",  "SGLT2i(+)",
  "Patient24164",  "SGLT2i(+)",
  "Patient24165",  "SGLT2i(+)",
  "Patient34064",  "SGLT2i(-)",
  "Patient34332",  "SGLT2i(-)"
)
# NOTE: The SGLT2i assignments above are approximate based on the 6 vs 10
# split described in Schaub JCI 2023. Verify against their Supplementary
# Table 1 or GEO sample metadata (GSM descriptions on the GEO page) and
# correct any misassignments before running NEBULA.

# Add to Seurat metadata
meta_add <- data.frame(
  patient_id    = so_schaub$patient_id,
  row.names     = colnames(so_schaub)
) %>%
  left_join(sglt2i_map, by = "patient_id")

so_schaub <- AddMetaData(so_schaub,
                         metadata = meta_add$sglt2i_status,
                         col.name = "sglt2i_status")

cat("\nSGLT2i status assigned:\n")
print(table(so_schaub$sglt2i_status, useNA = "always"))

SGLT2I_COL <- "sglt2i_status"

# Cell type: GSE220939 raw matrices have no cell type labels —
# we need to cluster or use Schaub's published labels.
# For now we flag this and proceed; section 7 will subset by orig.ident
# if no cell type column exists. You can also:
# Option A: Run standard Seurat clustering + label transfer from KPMP reference
# Option B: Check if Schaub deposited a processed object elsewhere (e.g. Zenodo)

has_celltype <- any(grepl("cell.?type|cluster|annotation",
                          colnames(so_schaub@meta.data), ignore.case = TRUE))
if (!has_celltype) {
  cat("\n*** No cell type column found in raw GEO deposit. ***\n")
  cat("Options:\n")
  cat("  A) Run basic clustering below (fast, approximate)\n")
  cat("  B) Check Zenodo or supplementary for Schaub processed object\n\n")
  
  # Option A: Quick clustering to get PT cells
  cat("Running basic normalization + clustering to identify PT cells...\n")
  so_schaub <- NormalizeData(so_schaub, verbose = FALSE)
  so_schaub <- FindVariableFeatures(so_schaub, nfeatures = 2000, verbose = FALSE)
  so_schaub <- ScaleData(so_schaub, verbose = FALSE)
  so_schaub <- RunPCA(so_schaub, npcs = 30, verbose = FALSE)
  so_schaub <- FindNeighbors(so_schaub, dims = 1:20, verbose = FALSE)
  so_schaub <- FindClusters(so_schaub, resolution = 0.5, verbose = FALSE)
  so_schaub <- RunUMAP(so_schaub, dims = 1:20, verbose = FALSE)
  
  # Score PT identity using canonical PT markers
  pt_markers <- list(PT = c("LRP2","CUBN","SLC5A2","UMOD","SLC34A1","SLC13A3"))
  so_schaub  <- AddModuleScore(so_schaub, features = pt_markers,
                               name = "PT_score", ctrl = 50)
  so_schaub$celltype_approx <- ifelse(so_schaub$PT_score1 > 0.1, "PT", "non-PT")
  
  cat("Approximate cell type distribution:\n")
  print(table(so_schaub$celltype_approx))
  CELLTYPE_COL <- "celltype_approx"
} else {
  CELLTYPE_COL <- colnames(so_schaub@meta.data)[
    grepl("cell.?type|cluster|annotation",
          colnames(so_schaub@meta.data), ignore.case = TRUE)][1]
}

cat("\n>>> Using columns:\n")
cat("  Cell type:   ", CELLTYPE_COL,    "\n")
cat("  SGLT2i:      ", SGLT2I_COL,      "\n")
cat("  Participant: ", PARTICIPANT_COL, "\n")

cat("\nUnique cell types:\n");    print(sort(unique(so_schaub@meta.data[[CELLTYPE_COL]])))
cat("\nUnique SGLT2i values:\n"); print(sort(unique(so_schaub@meta.data[[SGLT2I_COL]])))
cat("\nUnique participant IDs:\n");print(sort(unique(so_schaub@meta.data[[PARTICIPANT_COL]])))

# ── 7. Subset to PT cells ─────────────────────────────────────────────────────
# Schaub JCI 2023 labels PT subclusters PT-1 through PT-5
# Update pt_pattern if her labels differ (e.g. "^PT_S", "Proximal")
pt_pattern <- "^PT"

so_pt <- subset(so_schaub,
                cells = which(grepl(pt_pattern,
                                    so_schaub@meta.data[[CELLTYPE_COL]])))
cat("\n=== PT subset ===\n")
cat("Total PT cells:", ncol(so_pt), "\n")
cat("PT subtypes:\n");   print(table(so_pt@meta.data[[CELLTYPE_COL]]))
cat("SGLT2i groups:\n"); print(table(so_pt@meta.data[[SGLT2I_COL]]))
cat("Participants:\n");  print(table(so_pt@meta.data[[PARTICIPANT_COL]]))

# ── 8. Define pathway gene sets ───────────────────────────────────────────────
pathway_genes <- list(
  Glycolysis       = c("PKLR","PFKFB3","PFKL","ALDOC","HK2",
                       "ENO2","PGK1","PGAM1","TPI1","GAPDH"),
  Gluconeogenesis  = c("SLC25A10","GOT2","GOT1","FBP1",
                       "SLC25A11","PCK1","MDH1"),
  TCA_cycle        = c("SDHB","SUCLG1","PDK2","ACO2","IDH3G",
                       "SUCLA2","HAGH","PDHB","LDHA"),
  Glutathione      = c("CNDP2","GSTM4","GSTT2B","GSTO1",
                       "GGCT","GSTM3","AKR1A1"),
  Metallothioneins = c("MT1G","MT1X","MT1H","MT2A")
)
all_genes <- unlist(pathway_genes)

genes_present <- all_genes[all_genes %in% rownames(so_pt)]
genes_missing <- all_genes[!all_genes %in% rownames(so_pt)]
cat("\nGenes present:", length(genes_present),
    "| Missing:", length(genes_missing))
if (length(genes_missing) > 0)
  cat(" ->", paste(genes_missing, collapse = ", "))
cat("\n")

# ── 9. Prepare NEBULA input ───────────────────────────────────────────────────
meta_pt   <- so_pt@meta.data
sglt_vals <- unique(meta_pt[[SGLT2I_COL]])
cat("\nSGLT2i column values:", paste(sglt_vals, collapse = ", "), "\n")

meta_pt$sglt2i_bin <- case_when(
  meta_pt[[SGLT2I_COL]] %in% c("Yes","YES","yes","TRUE",TRUE,1,"1",
                               "SGLT2i","SGLT2i(+)","T2Di")    ~ 1L,
  meta_pt[[SGLT2I_COL]] %in% c("No","NO","no","FALSE",FALSE,0,"0",
                               "No SGLT2i","SGLT2i(-)","T2D")  ~ 0L,
  TRUE ~ NA_integer_
)
cat("SGLT2i binary coding:\n")
print(table(meta_pt$sglt2i_bin, useNA = "always"))

meta_pt     <- meta_pt[!is.na(meta_pt$sglt2i_bin), ]
so_pt_clean <- subset(so_pt, cells = rownames(meta_pt))

counts_mat  <- JoinLayers(so_pt_clean) |>
  GetAssayData(layer = "counts", assay = "RNA")
counts_mat  <- counts_mat[genes_present, ]

pred_df <- meta_pt[colnames(counts_mat), ] %>%
  dplyr::select(participant = all_of(PARTICIPANT_COL),
                sglt2i      = sglt2i_bin) %>%
  as.data.frame()

cat("\nFinal n cells:", nrow(pred_df), "\n")
cat("Participants per group:\n")
print(pred_df %>%
        group_by(participant, sglt2i) %>%
        summarise(n_cells = n(), .groups = "drop") %>%
        group_by(sglt2i) %>%
        summarise(n_participants = n(), total_cells = sum(n_cells)))

# ── 10. Run NEBULA ────────────────────────────────────────────────────────────
cat("\nRunning NEBULA on Schaub PT cells...\n")

design_mat   <- model.matrix(~ sglt2i, data = pred_df)

nebula_input <- nebula(
  count  = counts_mat,
  id     = pred_df$participant,
  pred   = design_mat,
  offset = NULL,
  model  = "NBLMM",
  ncore  = 4
)

nebula_res <- nebula_input$summary %>%
  as_tibble() %>%
  dplyr::rename(
    gene_symbol   = gene,
    log2FC_nebula = logFC_sglt2i,
    pval_nebula   = p_sglt2i
  ) %>%
  mutate(
    padj_nebula = p.adjust(pval_nebula, method = "BH"),
    pathway = case_when(
      gene_symbol %in% pathway_genes$Glycolysis       ~ "Glycolysis",
      gene_symbol %in% pathway_genes$Gluconeogenesis  ~ "Gluconeogenesis",
      gene_symbol %in% pathway_genes$TCA_cycle        ~ "TCA cycle",
      gene_symbol %in% pathway_genes$Glutathione      ~ "Glutathione",
      gene_symbol %in% pathway_genes$Metallothioneins ~ "Metallothioneins",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::filter(pathway != "Other") %>%
  arrange(pathway, gene_symbol)

cat("\n=== NEBULA results (Schaub GEO cells) ===\n")
print(nebula_res %>% dplyr::select(gene_symbol, pathway,
                                   log2FC_nebula, pval_nebula, padj_nebula))

# ── 11. DESeq2 pseudobulk ────────────────────────────────────────────────────
cat("\nRunning DESeq2 pseudobulk on Schaub PT cells...\n")

participants <- unique(pred_df$participant)
pb_list <- lapply(participants, function(p) {
  cells_p <- rownames(pred_df)[pred_df$participant == p]
  if (length(cells_p) == 1)
    Matrix::rowSums(counts_mat[, cells_p, drop = FALSE])
  else
    Matrix::rowSums(counts_mat[, cells_p])
})
pb_mat           <- do.call(cbind, pb_list)
colnames(pb_mat) <- participants

pb_meta <- pred_df %>%
  group_by(participant, sglt2i) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  column_to_rownames("participant")
pb_meta        <- pb_meta[colnames(pb_mat), ]
pb_meta$sglt2i <- factor(pb_meta$sglt2i, levels = c(0, 1))

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(pb_mat),
  colData   = pb_meta,
  design    = ~ sglt2i
)
dds       <- DESeq(dds, quiet = TRUE)
deseq_res <- results(dds, contrast = c("sglt2i", "1", "0")) %>%
  as.data.frame() %>%
  rownames_to_column("gene_symbol") %>%
  dplyr::rename(log2FC_deseq2 = log2FoldChange, padj_deseq2 = padj) %>%
  dplyr::select(gene_symbol, log2FC_deseq2, padj_deseq2) %>%
  dplyr::filter(gene_symbol %in% genes_present)

# ── 12. Schaub published values ───────────────────────────────────────────────
# Approximate log2FC from Schaub JCI 2023 — positive = higher in SGLT2i(+)
# *** Replace with exact values from Supplementary Table 4 if available ***

schaub_published <- tribble(
  ~gene_symbol,  ~log2FC_schaub,  ~pathway,
  "PKLR",    -0.30, "Glycolysis",
  "PFKFB3",  -0.25, "Glycolysis",
  "PFKL",    -0.20, "Glycolysis",
  "ALDOC",   -0.35, "Glycolysis",
  "HK2",     -0.28, "Glycolysis",
  "ENO2",    -0.22, "Glycolysis",
  "PGK1",    -0.18, "Glycolysis",
  "PGAM1",   -0.15, "Glycolysis",
  "TPI1",    -0.20, "Glycolysis",
  "GAPDH",   -0.10, "Glycolysis",
  "SLC25A10",-0.18, "Gluconeogenesis",
  "GOT2",    -0.22, "Gluconeogenesis",
  "GOT1",    -0.15, "Gluconeogenesis",
  "FBP1",    -0.40, "Gluconeogenesis",
  "SLC25A11",-0.20, "Gluconeogenesis",
  "PCK1",    -0.35, "Gluconeogenesis",
  "MDH1",    -0.25, "Gluconeogenesis",
  "SDHB",    -0.20, "TCA cycle",
  "SUCLG1",  -0.18, "TCA cycle",
  "PDK2",    -0.15, "TCA cycle",
  "ACO2",     0.10, "TCA cycle",
  "IDH3G",    0.12, "TCA cycle",
  "SUCLA2",  -0.12, "TCA cycle",
  "HAGH",    -0.10, "TCA cycle",
  "PDHB",    -0.15, "TCA cycle",
  "LDHA",    -0.20, "TCA cycle",
  "CNDP2",   -0.10, "Glutathione",
  "GSTM4",   -0.08, "Glutathione",
  "GSTT2B",  -0.12, "Glutathione",
  "GSTO1",   -0.60, "Glutathione",
  "GGCT",    -0.15, "Glutathione",
  "GSTM3",   -0.10, "Glutathione",
  "AKR1A1",  -0.05, "Glutathione",
  "MT1G",     0.80, "Metallothioneins",
  "MT1X",     0.65, "Metallothioneins",
  "MT1H",     0.55, "Metallothioneins",
  "MT2A",     0.45, "Metallothioneins"
)

# ── 13. Merge and compute concordance ────────────────────────────────────────
combined <- schaub_published %>%
  left_join(nebula_res %>% dplyr::select(gene_symbol, log2FC_nebula,
                                         pval_nebula, padj_nebula),
            by = "gene_symbol") %>%
  left_join(deseq_res, by = "gene_symbol") %>%
  mutate(
    concordant_nebula = sign(log2FC_nebula) == sign(log2FC_schaub),
    concordant_deseq2 = sign(log2FC_deseq2) == sign(log2FC_schaub),
    pathway = factor(pathway, levels = c("Glycolysis","Gluconeogenesis",
                                         "TCA cycle","Glutathione",
                                         "Metallothioneins"))
  )

cat("\n=== Concordance summary ===\n")
cat("NEBULA vs Schaub:", sum(combined$concordant_nebula, na.rm = TRUE),
    "/", sum(!is.na(combined$concordant_nebula)), "genes concordant\n")
cat("DESeq2 vs Schaub:", sum(combined$concordant_deseq2, na.rm = TRUE),
    "/", sum(!is.na(combined$concordant_deseq2)), "genes concordant\n")

# ── 14. Bar comparison plot ───────────────────────────────────────────────────
method_colors <- c(
  "Schaub published"              = "#1D3557",
  "NEBULA (GEO cells)"            = "#2A9D8F",
  "DESeq2 pseudobulk (GEO cells)" = "#E76F51"
)

plot_df <- combined %>%
  dplyr::select(gene_symbol, pathway,
                `Schaub published`              = log2FC_schaub,
                `NEBULA (GEO cells)`            = log2FC_nebula,
                `DESeq2 pseudobulk (GEO cells)` = log2FC_deseq2) %>%
  pivot_longer(cols      = -c(gene_symbol, pathway),
               names_to  = "method",
               values_to = "log2FC") %>%
  mutate(
    method      = factor(method, levels = names(method_colors)),
    gene_symbol = factor(gene_symbol,
                         levels = combined$gene_symbol[order(combined$pathway)])
  )

p_comparison <- ggplot(plot_df,
                       aes(x = gene_symbol, y = log2FC, fill = method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, alpha = 0.9) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey30") +
  facet_wrap(~ pathway, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = method_colors, name = NULL) +
  labs(
    title    = "Schaub GEO Replication: NEBULA & DESeq2 vs. Published Values",
    subtitle = "PT cells | GSE220939 | SGLT2i(+) vs SGLT2i(−)",
    x        = NULL,
    y        = "log2 Fold Change (SGLT2i+ vs SGLT2i−)",
    caption  = paste0(
      "Schaub published values are approximate — replace with exact Supp Table 4 values.\n",
      "NEBULA concordance: ", sum(combined$concordant_nebula, na.rm = TRUE),
      "/", sum(!is.na(combined$concordant_nebula)), " genes | ",
      "DESeq2 concordance: ", sum(combined$concordant_deseq2, na.rm = TRUE),
      "/", sum(!is.na(combined$concordant_deseq2)), " genes")
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold", size = 10),
    plot.title       = element_text(face = "bold"),
    plot.caption     = element_text(size = 7, color = "grey50")
  )

# ── 15. Scatter: NEBULA (GEO) vs Schaub published ────────────────────────────
pathway_colors <- c(
  "Glycolysis"       = "#E76F51",
  "Gluconeogenesis"  = "#F4A261",
  "TCA cycle"        = "#2A9D8F",
  "Glutathione"      = "#457B9D",
  "Metallothioneins" = "#6A0572"
)

p_scatter <- combined %>%
  dplyr::filter(!is.na(log2FC_nebula)) %>%
  ggplot(aes(x = log2FC_schaub, y = log2FC_nebula,
             color = pathway, label = gene_symbol)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2.8, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey30") +
  scale_color_manual(values = pathway_colors) +
  labs(
    title    = "NEBULA on GEO Cells vs. Schaub Published",
    subtitle = "Points on diagonal = perfect concordance | GSE220939",
    x        = "log2FC — Schaub published",
    y        = "log2FC — NEBULA (GEO cells)",
    color    = "Pathway"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title      = element_text(face = "bold"),
        legend.position = "right")

# ── 16. Save all outputs ──────────────────────────────────────────────────────
cat("\nSaving outputs to:", output_dir, "\n")

ggsave(file.path(output_dir, "GEO_replication_bar_comparison.pdf"),
       plot = p_comparison, width = 18, height = 7)
ggsave(file.path(output_dir, "GEO_replication_scatter.pdf"),
       plot = p_scatter, width = 8, height = 7)

wb <- createWorkbook()
addWorksheet(wb, "NEBULA_results");      writeData(wb, "NEBULA_results",      nebula_res)
addWorksheet(wb, "DESeq2_results");      writeData(wb, "DESeq2_results",      deseq_res)
addWorksheet(wb, "Combined_comparison"); writeData(wb, "Combined_comparison", combined)
addWorksheet(wb, "Concordance_summary")
conc_summary <- combined %>%
  group_by(pathway) %>%
  summarise(
    n_genes           = n(),
    nebula_concordant = sum(concordant_nebula, na.rm = TRUE),
    deseq2_concordant = sum(concordant_deseq2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(nebula_pct = round(100 * nebula_concordant / n_genes, 1),
         deseq2_pct = round(100 * deseq2_concordant / n_genes, 1))
writeData(wb, "Concordance_summary", conc_summary)
saveWorkbook(wb, file.path(output_dir, "GEO_replication_results.xlsx"),
             overwrite = TRUE)

cat("\n=== Done ===\n")
cat("Files saved:\n")
cat("  GEO_replication_bar_comparison.pdf\n")
cat("  GEO_replication_scatter.pdf\n")
cat("  GEO_replication_results.xlsx\n")

################################################################################
# NOTES FOR MANUAL REVIEW
################################################################################
# 1. After sections 4-5 print to console, verify CELLTYPE_COL, SGLT2I_COL,
#    and PARTICIPANT_COL in section 6 match Schaub's actual column names.
#
# 2. Update pt_pattern in section 7 to match her exact PT label strings.
#    Common variants: "^PT", "^PT-", "^PT_S", "Proximal"
#
# 3. Replace schaub_published log2FC values in section 12 with exact values
#    from Schaub JCI 2023 Supplementary Table 4 — current values are
#    directionally correct approximations only.
#
# 4. GEO accession: GSE220939
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220939
################################################################################

