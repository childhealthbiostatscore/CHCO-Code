# merge PB90 + ATTEMPT scRNA

library(Seurat)
library(jsonlite)
library(aws.s3)
library(dplyr)
library(purrr)
user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
} else if (user == "yejichoi") { # hyak version
  root_path <- "/mmfs1/gscratch/togo/yejichoi"
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
} else {
  stop("Unknown user: please specify root path for this user.")
}

## Create an S3 client
Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

s3read_using_region <- function(FUN, ..., object, bucket, region = NULL, opts = NULL, filename = NULL) {
  if (missing(bucket)) {
    bucket <- get_bucketname(object)
  }
  object <- get_objectkey(object)
  
  tmp <- if (is.character(filename)) {
    file.path(tempdir(TRUE), filename)
  } else {
    tempfile(fileext = paste0(".", tools::file_ext(object)))
  }
  
  on.exit(unlink(tmp))
  
  # Add region to opts if provided
  if (!is.null(region)) {
    if (is.null(opts)) {
      opts <- list(region = region)
    } else {
      opts$region <- region
    }
  }
  
  if (is.null(opts)) {
    r <- save_object(bucket = bucket, object = object, file = tmp)
  } else {
    r <- do.call("save_object", c(list(bucket = bucket, 
                                       object = object, 
                                       file = tmp), opts))
  }
  
  return(FUN(tmp, ...))
}

# pull data
pb90 <- s3readRDS("Kidney transcriptomics/Single cell RNA seq/PB_90_RPCAFix_Metadata_LR_RM_kpmpV1labelled.rds", "scrna", region = "")
attempt <- s3readRDS("Kidney transcriptomics/Single cell RNA seq/PB_attempt_harmony_rpca_RM_kpmpV1labelled_Sept2024.RDS", "scrna", region = "")

# attempt_dat <- s3readRDS("attempt_clinical_data.RDS", "harmonized.dataset", region = "")
# 
# attempt_meta <- attempt@meta.data %>%
#   dplyr::mutate(subject_id = Subject.ID,
#                 visit = case_when(Visit == "baseline" ~ "PRE", 
#                                   Visit == "post_treatment" ~ "POST"))
# attempt_meta <- left_join(attempt_meta, attempt_dat, by = c("subject_id", "visit"))
# rownames(attempt_meta) <- attempt_meta$barcode
# attempt <- AddMetaData(attempt, attempt_meta)
# View(attempt@meta.data)
# add/modify variables in metadata for easy integration
colnames(attempt@meta.data)[colnames(attempt@meta.data) == "Study"] <- "cohort"
colnames(attempt@meta.data)[colnames(attempt@meta.data) == "Subject.ID"] <- "record_id"
colnames(attempt@meta.data)[colnames(attempt@meta.data) == "Co.Enrolled"] <- "co.enroll.id"
colnames(attempt@meta.data)[colnames(attempt@meta.data) == "Sex"] <- "sex"
colnames(attempt@meta.data)[colnames(attempt@meta.data) == "Visit"] <- "visit"
colnames(attempt@meta.data)[colnames(attempt@meta.data) == "Dx"] <- "group"
colnames(attempt@meta.data)[colnames(attempt@meta.data) == "celltype"] <- "celltype_harmony"
colnames(attempt@meta.data)[colnames(attempt@meta.data) == "Cryostor_ID"] <- "cryostor_id"
attempt$visit <- ifelse(attempt$visit == "BL", "baseline", "post_treatment")
attempt$X <- NULL
attempt$Site <- NULL
attempt$source = "attempt"

# attempt trt assignment
attempt_grp <- s3read_using_region(bucket = "attempt", object = "Clinical Data/ATTEMPT_unblinding.csv", FUN = read.csv, region = "") %>%
  dplyr::select(record_id = subject_id, treatment)
attempt@meta.data <- attempt@meta.data %>% left_join(attempt_grp, by = "record_id")
rownames(attempt@meta.data) <- attempt$barcode
pb90$group <- gsub("_", " ", pb90$group)
pb90$source = "pb90"

pb90_attempt <- merge(pb90, y = attempt, 
                      add.cell.ids = c("pb90", "attempt"),
                      project = "pb90_attempt")

# add meta data
harm_dat <- s3read_using_region(FUN=read.csv, na = "", object = "harmonized_dataset.csv", bucket = "harmonized.dataset", region = "")
harm_dat <- harm_dat %>%
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(.x, na.rm = TRUE))),
                   .by = c(record_id, visit))
pb90_attempt_meta <- pb90_attempt@meta.data %>%
  dplyr::select(-c(cryostor_id, kit_id, group, sex, age, sglt2i_ever, hba1c, eGFR_CKD_epi))
pb90_attempt_meta <- left_join(pb90_attempt_meta, harm_dat, by = c("record_id", "visit"))
rownames(pb90_attempt_meta) <- pb90_attempt_meta$barcode
pb90_attempt <- AddMetaData(pb90_attempt, pb90_attempt_meta)

# normalize/scale
pb90_attempt <- NormalizeData(pb90_attempt)
pb90_attempt <- FindVariableFeatures(pb90_attempt)
pb90_attempt <- ScaleData(pb90_attempt)
pb90_attempt <- RunPCA(pb90_attempt)

pb90_attempt <- FindNeighbors(pb90_attempt, dims = 1:30, reduction = "pca")
pb90_attempt <- FindClusters(pb90_attempt, resolution = 1.5, cluster.name = "unintegrated_clusters")

pb90_attempt <- RunUMAP(pb90_attempt, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# Data Integration using Harmony
pb90_attempt_integrated <- IntegrateLayers(object = pb90_attempt, method = HarmonyIntegration, 
                                           orig.reduction = "pca", new.reduction = "harmony", verbose = T)

pb90_attempt_integrated <- FindNeighbors(pb90_attempt_integrated, reduction = "harmony", dims = 1:30)
pb90_attempt_integrated <- FindClusters(pb90_attempt_integrated, resolution = 1.5, cluster.name = "harmony_clusters")

pb90_attempt_integrated <- RunUMAP(pb90_attempt_integrated, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

p1 <- DimPlot(pb90_attempt, reduction = "umap.unintegrated", group.by = c("source"))
p2 <- DimPlot(pb90_attempt_integrated, reduction = "umap.harmony", group.by = c("source"))
p1 | p2

# # Join all layers of the same type
pb90_attempt_integrated <- SeuratObject::JoinLayers(pb90_attempt_integrated)
# Remove all scale.data layers
pb90_attempt_integrated[["RNA"]]$scale.data.1 <- NULL
pb90_attempt_integrated[["RNA"]]$scale.data.2 <- NULL
pb90_attempt_integrated[["RNA"]]$scale.data <- NULL

# Re-scale with the joined data
pb90_attempt_integrated <- ScaleData(pb90_attempt_integrated)

s3saveRDS(pb90_attempt_integrated, "data_clean/pb90_attempt_integrated.rds", "scrna", region = "", multipart = T)

# Process data for downstream analysis
ncol(pb90_attempt_integrated) # Number of cells before filtering
nrow(pb90_attempt_integrated) # Number of genes before filtering

# Filter out rare genes expressed in less than "gene_pct" of cells
expr_matrix <- as.matrix(GetAssayData(pb90_attempt_integrated, assay = "RNA", layer = "counts"))

# Calculate the proportion of cells expressing each gene
num_cells_per_gene <- rowSums(expr_matrix > 0)  # Count nonzero cells per gene
total_cells <- ncol(expr_matrix)  # Total number of cells
gene_proportion <- num_cells_per_gene / total_cells  # Fraction of cells expressing each gene

# Keep genes expressed in at least "gene_pct" of cells
gene_pct = 0.05 # 5%, change based on biopsy type & data
genes_to_keep <- names(gene_proportion[gene_proportion >= gene_pct])
pb90_attempt_integrated <- subset(pb90_attempt_integrated, features = genes_to_keep)

#Check cells & genes remaining
ncol(pb90_attempt_integrated) # Number of cells remaining after filtering
nrow(pb90_attempt_integrated) # Number of genes remaining after filtering

#Check the number of Mitochondrial genes to start
sum(grepl("^MT-", rownames(pb90_attempt_integrated))) 

# Identify mitochondrial genes (human: start with "MT-")
mito_genes <- grep("^MT-", rownames(pb90_attempt_integrated), value = TRUE)
pb90_attempt_integrated <- subset(pb90_attempt_integrated, features = setdiff(rownames(pb90_attempt_integrated), mito_genes))

#Check the number of Mitochondrial genes after filtering to ensure filtering step was successful
sum(grepl("^MT-", rownames(pb90_attempt_integrated))) #Should be 0

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

pb90_attempt_integrated <- subset(pb90_attempt_integrated, features = setdiff(rownames(pb90_attempt_integrated), ribo_genes))

s3saveRDS(pb90_attempt_integrated, "data_clean/pb90_attempt_integrated_processed.rds", "scrna", region = "", multipart = T)
