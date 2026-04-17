library(dplyr)
library(purrr)
library(furrr)
library(aws.s3)
library(jsonlite)
library(ggplot2)
library(ggforce)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(msigdbr)
library(ggpubr)
library(table1)
library(nlme)
library(reshape2)
library(Seurat)
library(nebula)
library(doParallel)

user <- Sys.info()[["user"]]

if (user == "choiyej") {
  # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON(
    "/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json"
  )
} else if (user == "rameshsh") {
  # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/rameshsh/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") {
  # hyak version
  root_path <- "/mmfs1/gscratch/togo/yejichoi/"
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/pylell/keys.json")
} else if (user == "kristenmiller") {
  root_path <- "/Users/kristenmiller/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive"
  git_path <- "/Users/kristenmiller/Documents/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/Users/kristenmiller/keys.json")
} else {
  stop("Unknown user: please specify root path for this user.")
}

source(file.path(git_path, "Renal HEIRitage/RH_RH2_IMPROVE_functions.R"))
source(file.path(git_path, "Renal HEIRitage/RH_RH2_IMPROVE_scRNA_functions.R"))

entrez_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
  columns = c("SYMBOL", "GENENAME", "ALIAS"),
  keytype = "ENTREZID"
)

s3write_using_region <- function(
  FUN,
  ...,
  object,
  bucket,
  region = NULL,
  opts = NULL,
  filename = NULL
) {
  if (missing(bucket)) {
    bucket <- get_bucketname(object)
  }
  object <- get_objectkey(object)

  tmp <- if (is.character(filename)) {
    file.path(tempdir(TRUE), filename)
  } else {
    # if object has an extension, keep it; otherwise make a generic tmp
    ext <- tools::file_ext(object)
    if (nzchar(ext)) tempfile(fileext = paste0(".", ext)) else tempfile()
  }

  on.exit(unlink(tmp), add = TRUE)

  # Add region to opts if provided
  if (!is.null(region)) {
    if (is.null(opts)) {
      opts <- list(region = region)
    } else {
      opts$region <- region
    }
  }

  FUN(tmp, ...)

  if (is.null(opts)) {
    r <- put_object(file = tmp, bucket = bucket, object = object)
  } else {
    r <- do.call(
      "put_object",
      c(list(file = tmp, bucket = bucket, object = object), opts)
    )
  }

  return(invisible(r))
}

# read in gene lists from Christine's excel file:
gene_lists <- openxlsx::read.xlsx(
  "./Limonte/Autophagy_gene_protein_lists_CL.xlsx",
  sheet = "Gene lists",
  startRow = 3,
  colNames = TRUE
)

# read in full pb90 seurat object:
pb90_so <- s3readRDS("data_clean/pb90_processed.rds", "scrna", region = "")

# remove obese controls:
pb90_sub <- subset(
  pb90_so,
  group %in% c("Lean Control", "Type 1 Diabetes", "Type 2 Diabetes")
)
# only use baseline samples:
pb90_sub <- subset(pb90_sub, visit %in% c("baseline"))

# save two objects, one for all cells, one for PT only
pb90_allcells <- pb90_sub
# only use proximal tubular cells:
pb90_ptcells <- subset(
  pb90_sub,
  KPMP_celltype %in% c("aPT", "PT-S1/S2", "PT-S3")
)

#create gene lists from Christine's files
overall_autophagy_genes <- c(
  gene_lists$`Overall.autophagy-related.genes.(n=589)`[
    !is.na(gene_lists$`Overall.autophagy-related.genes.(n=589)`)
  ]
)
autophagy_induction_genes <- c(
  gene_lists$`Autophagy.induction.genes.(n=20)`[
    !is.na(gene_lists$`Autophagy.induction.genes.(n=20)`)
  ]
)
lysosomal_biogenesis_genes <- c(
  gene_lists$`Lysosomal.biogenesis.genes.(n=26)`[
    !is.na(gene_lists$`Lysosomal.biogenesis.genes.(n=26)`)
  ]
)
tfeb_genes <- c(
  gene_lists$`TFEB.genes.(n=13)`[!is.na(gene_lists$`TFEB.genes.(n=13)`)]
)
genes_list <- c(
  overall_autophagy_genes,
  autophagy_induction_genes,
  lysosomal_biogenesis_genes,
  tfeb_genes
)

# nebula for gene level analysis:
#8.2.1 Prepare Data for Analysis

#### (1): Diabetes vs. Healthy Controls; PT cells ####
# Load counts data from the object and round to nearest integer
# use counts_hvg naming to match pipeline
seurat_object_hvg <- pb90_ptcells
counts_hvg <- round(GetAssayData(seurat_object_hvg, layer = "counts"))

# Create a list of genes from the row names of the counts matrix
# Make sure the exposure/independent/x variable (group variable) is a factor
seurat_object_hvg$group <- as.factor(ifelse(
  seurat_object_hvg$group == "Lean Control",
  0,
  1
))

source("/Users/tim/Downloads/nebula_pipeline.R")
# Nebula all genes together
seuratdata = as.SingleCellExperiment(seurat_object_hvg)
seuratdata <- scToNeb(
  obj = seuratdata,
  assay = "RNA",
  id = "record_id",
  pred = "group",
  offset = "nCount_RNA"
)
seuratdata$count = round(seuratdata$count)
seuratdata$pred = model.matrix(~group, seuratdata$pred)
re = nebula(
  seuratdata$count,
  seuratdata$id,
  pred = seuratdata$pred,
  ncore = 8
)
diabetes_ptcell_nebula_results_all = re$summary[
  match(full_results$gene, re$summary$gene),
]

# save results:
diabetes_ptcell_nebula_results <- full_results
save(
  diabetes_ptcell_nebula_results,
  diabetes_ptcell_nebula_results_all,
  file = "/Users/tim/Downloads/diabetes_ptcell_nebula_results.RData"
)
rm(diabetes_ptcell_nebula_results, diabetes_ptcell_nebula_results_all)
# ---- Volcano Plot ----
# Create the volcano plot using ggplot
diabetes_ptcell_volcano_plot <- ggplot(
  full_results,
  aes(x = logFC_group1, y = PValue10, color = color)
) +
  geom_point(alpha = 0.7) + # Add points with transparency
  scale_color_identity() + # Use color directly from 'color' column
  theme_minimal() + # Apply minimalistic theme
  labs(
    title = "Diabetes vs. Healthy Controls", # Main title
    subtitle = "PT Cells Only", # Subtitle
    x = "logFC", # X-axis label
    y = "-log10(P-Value)", # Y-axis label
    color = "LogFC Direction",
    caption = paste0(
      "FDR < 0.05, Genes = ",
      Genes,
      ", Nuclei = ",
      Nuclei,
      ", Non-Convergence Rate: ",
      Nonconvergence_Rate,
      ", Genes Filtered out for Low Expression: ",
      low_exp
    ) # Plot caption
  ) +
  theme(
    plot.title = element_text(hjust = 0),
    axis.text.x = element_text(angle = 0, hjust = 1)
  ) +
  xlim(min, max) + # Set x-axis limits
  # Add gene labels for significant points
  geom_text(
    data = significant_df,
    aes(label = gene),
    vjust = 1,
    hjust = 1,
    size = 3,
    check_overlap = TRUE,
    color = "black"
  )

# Display the volcano plot
diabetes_ptcell_volcano_plot
ggsave(
  "/Users/tim/Downloads/diabetes_ptcell_volcano_plot.png",
  plot = diabetes_ptcell_volcano_plot,
  dpi = 300
)

#### (2): SGLT2i vs. no; PT cells ####
# Load counts data from the object and round to nearest integer
# use counts_hvg to match pipeline
pb90_ptcells_t2d <- subset(
  pb90_ptcells,
  pb90_ptcells$group == "Type 2 Diabetes"
)
seurat_object_hvg <- pb90_ptcells_t2d
counts_hvg <- round(GetAssayData(seurat_object_hvg, layer = "counts"))

# Create a list of genes from the row names of the counts matrix
genes_list <- c(
  overall_autophagy_genes,
  autophagy_induction_genes,
  lysosomal_biogenesis_genes,
  tfeb_genes
)

# Make sure the exposure/independent/x variable (group variable) is a factor
seurat_object_hvg$group <- as.factor(ifelse(
  seurat_object_hvg$sglt2i_ever == "No",
  0,
  1
))

source("/Users/tim/Downloads/nebula_pipeline.R")
sglt2i_ptcell_nebula_results <- full_results
save(
  sglt2i_ptcell_nebula_results,
  file = "/Users/tim/Downloads/sglt2i_ptcell_nebula_results.RData"
)
rm(sglt2i_ptcell_nebula_results)
# ---- Volcano Plot ----
# Create the volcano plot using ggplot
sglt2i_ptcell_volcano_plot <- ggplot(
  full_results,
  aes(x = logFC_group1, y = PValue10, color = color)
) +
  geom_point(alpha = 0.7) + # Add points with transparency
  scale_color_identity() + # Use color directly from 'color' column
  theme_minimal() + # Apply minimalistic theme
  labs(
    title = "SGLT2i vs. No SGLT2i", # Main title
    subtitle = "PT Cells Only", # Subtitle
    x = "logFC", # X-axis label
    y = "-log10(P-Value)", # Y-axis label
    color = "LogFC Direction",
    caption = paste0(
      "FDR < 0.05, Genes = ",
      Genes,
      ", Nuclei = ",
      Nuclei,
      ", Non-Convergence Rate: ",
      Nonconvergence_Rate,
      ", Genes Filtered out for Low Expression: ",
      low_exp
    ) # Plot caption
  ) +
  theme(
    plot.title = element_text(hjust = 0),
    axis.text.x = element_text(angle = 0, hjust = 1)
  ) +
  xlim(min, max) + # Set x-axis limits
  # Add gene labels for significant points
  geom_text(
    data = significant_df,
    aes(label = gene),
    vjust = 1,
    hjust = 1,
    size = 3,
    check_overlap = TRUE,
    color = "black"
  )

# Display the volcano plot
sglt2i_ptcell_volcano_plot
ggsave(
  "/Users/tim/Downloads/sglt2i_ptcell_volcano_plot.png",
  plot = sglt2i_ptcell_volcano_plot,
  dpi = 300
)

#### (3): Diabetes vs. Healthy Controls; all cells ####
# Load counts data from the object and round to nearest integer
# use counts_hvg to match pipeline
seurat_object_hvg <- pb90_allcells
counts_hvg <- round(GetAssayData(seurat_object_hvg, layer = "counts"))

# Create a list of genes from the row names of the counts matrix
genes_list <- c(
  overall_autophagy_genes,
  autophagy_induction_genes,
  lysosomal_biogenesis_genes,
  tfeb_genes
)

# Make sure the exposure/independent/x variable (group variable) is a factor
seurat_object_hvg$group <- as.factor(ifelse(
  seurat_object_hvg$group == "Lean Control",
  0,
  1
))

source("/Users/tim/Downloads/nebula_pipeline.R")
diabetes_allcell_nebula_results <- full_results
save(
  diabetes_allcell_nebula_results,
  file = "/Users/tim/Downloads/diabetes_allcell_nebula_results.RData"
)
rm(diabetes_allcell_nebula_results)
# ---- Volcano Plot ----
# Create the volcano plot using ggplot
diabetes_allcell_volcano_plot <- ggplot(
  full_results,
  aes(x = logFC_group1, y = PValue10, color = color)
) +
  geom_point(alpha = 0.7) + # Add points with transparency
  scale_color_identity() + # Use color directly from 'color' column
  theme_minimal() + # Apply minimalistic theme
  labs(
    title = "Diabetes vs. Healthy Controls", # Main title
    subtitle = "All cells", # Subtitle
    x = "logFC", # X-axis label
    y = "-log10(P-Value)", # Y-axis label
    color = "LogFC Direction",
    caption = paste0(
      "FDR < 0.05, Genes = ",
      Genes,
      ", Nuclei = ",
      Nuclei,
      ", Non-Convergence Rate: ",
      Nonconvergence_Rate,
      ", Genes Filtered out for Low Expression: ",
      low_exp
    ) # Plot caption
  ) +
  theme(
    plot.title = element_text(hjust = 0),
    axis.text.x = element_text(angle = 0, hjust = 1)
  ) +
  xlim(min, max) + # Set x-axis limits
  # Add gene labels for significant points
  geom_text(
    data = significant_df,
    aes(label = gene),
    vjust = 1,
    hjust = 1,
    size = 3,
    check_overlap = TRUE,
    color = "black"
  )

# Display the volcano plot
diabetes_allcell_volcano_plot
ggsave(
  "/Users/tim/Downloads/diabetes_allcell_volcano_plot.png",
  plot = diabetes_allcell_volcano_plot,
  dpi = 300
)

#### (2): SGLT2i vs. no; all cells ####
# Load counts data from the object and round to nearest integer
# use counts_hvg to match pipeline
pb90_allcells_t2d <- subset(
  pb90_allcells,
  pb90_allcells$group == "Type 2 Diabetes"
)
seurat_object_hvg <- pb90_allcells_t2d
counts_hvg <- round(GetAssayData(seurat_object_hvg, layer = "counts"))

# Create a list of genes from the row names of the counts matrix
genes_list <- c(
  overall_autophagy_genes,
  autophagy_induction_genes,
  lysosomal_biogenesis_genes,
  tfeb_genes
)

# Make sure the exposure/independent/x variable (group variable) is a factor
seurat_object_hvg$group <- as.factor(ifelse(
  seurat_object_hvg$sglt2i_ever == "No",
  0,
  1
))

source("/Users/tim/Downloads/nebula_pipeline.R")
sglt2i_allcell_nebula_results <- full_results
save(
  sglt2i_allcell_nebula_results,
  file = "/Users/tim/Downloads/sglt2i_allcell_nebula_results.RData"
)
rm(sglt2i_allcell_nebula_results)
# ---- Volcano Plot ----
# Create the volcano plot using ggplot
sglt2i_allcell_volcano_plot <- ggplot(
  full_results,
  aes(x = logFC_group1, y = PValue10, color = color)
) +
  geom_point(alpha = 0.7) + # Add points with transparency
  scale_color_identity() + # Use color directly from 'color' column
  theme_minimal() + # Apply minimalistic theme
  labs(
    title = "SGLT2i vs. No SGLT2i", # Main title
    subtitle = "All cells", # Subtitle
    x = "logFC", # X-axis label
    y = "-log10(P-Value)", # Y-axis label
    color = "LogFC Direction",
    caption = paste0(
      "FDR < 0.05, Genes = ",
      Genes,
      ", Nuclei = ",
      Nuclei,
      ", Non-Convergence Rate: ",
      Nonconvergence_Rate,
      ", Genes Filtered out for Low Expression: ",
      low_exp
    ) # Plot caption
  ) +
  theme(
    plot.title = element_text(hjust = 0),
    axis.text.x = element_text(angle = 0, hjust = 1)
  ) +
  xlim(min, max) + # Set x-axis limits
  # Add gene labels for significant points
  geom_text(
    data = significant_df,
    aes(label = gene),
    vjust = 1,
    hjust = 1,
    size = 3,
    check_overlap = TRUE,
    color = "black"
  )

# Display the volcano plot
sglt2i_allcell_volcano_plot
ggsave(
  "/Users/tim/Downloads/sglt2i_allcell_volcano_plot.png",
  plot = sglt2i_allcell_volcano_plot,
  dpi = 300
)
