# Rename rownames in Seurat object from gene IDs to gene names
# Assumes: sc_org is your Seurat object
# gene_ids <- rownames(sc_org)
# new_gene_names <- sc_org@misc$gene_annotations$name

library(Seurat)

# Function to rename Seurat object rownames
rename_seurat_genes <- function(seurat_obj, new_gene_names) {
  
  # Get current gene IDs
  gene_ids <- rownames(seurat_obj)
  
  # Verify lengths match
  if (length(gene_ids) != length(new_gene_names)) {
    stop("Length of gene_ids and new_gene_names must match")
  }
  
  # Handle potential duplicated gene names by making them unique
  new_gene_names_unique <- make.unique(new_gene_names)
  
  # Create a copy of the original object
  new_seurat <- seurat_obj
  
  # Get assay names
  assay_names <- names(seurat_obj@assays)
  
  # Process each assay
  for (assay_name in assay_names) {
    # Access assay using $ operator
    if (assay_name == "RNA") {
      assay <- seurat_obj@assays$RNA
    } else if (assay_name == "ADT") {
      assay <- seurat_obj@assays$ADT
    } else if (assay_name == "HTO") {
      assay <- seurat_obj@assays$HTO
    } else {
      # For other assays, use slot() function as backup
      assay <- slot(seurat_obj@assays, assay_name)
    }
    
    # Rename rownames in counts matrix
    if (!is.null(assay@counts) && nrow(assay@counts) > 0) {
      rownames(assay@counts) <- new_gene_names_unique
    }
    
    # Rename rownames in data matrix
    if (!is.null(assay@data) && nrow(assay@data) > 0) {
      rownames(assay@data) <- new_gene_names_unique
    }
    
    # Rename rownames in scale.data matrix
    if (!is.null(assay@scale.data) && nrow(assay@scale.data) > 0) {
      # Only rename the genes that are present in scale.data
      scale_data_genes <- rownames(assay@scale.data)
      gene_mapping <- setNames(new_gene_names_unique, gene_ids)
      new_scale_names <- gene_mapping[scale_data_genes]
      rownames(assay@scale.data) <- new_scale_names
    }
    
    # Update meta.features - this is crucial!
    if (!is.null(assay@meta.features) && nrow(assay@meta.features) > 0) {
      rownames(assay@meta.features) <- new_gene_names_unique
    }
    
    # Update variable features if they exist
    if (length(assay@var.features) > 0) {
      # Map variable features from old names to new names
      gene_mapping <- setNames(new_gene_names_unique, gene_ids)
      var_features_new <- gene_mapping[assay@var.features]
      var_features_new <- var_features_new[!is.na(var_features_new)]
      assay@var.features <- as.character(var_features_new)
    }
    
    # Put the updated assay back
    if (assay_name == "RNA") {
      new_seurat@assays$RNA <- assay
    } else if (assay_name == "ADT") {
      new_seurat@assays$ADT <- assay
    } else if (assay_name == "HTO") {
      new_seurat@assays$HTO <- assay
    } else {
      slot(new_seurat@assays, assay_name) <- assay
    }
  }
  
  # Restore dimensional reductions if they exist
  if (length(seurat_obj@reductions) > 0) {
    for (reduction_name in names(seurat_obj@reductions)) {
      reduction <- seurat_obj@reductions[[reduction_name]]
      
      # Update feature loadings if they exist and contain gene names
      if (!is.null(reduction@feature.loadings) && nrow(reduction@feature.loadings) > 0) {
        old_loading_names <- rownames(reduction@feature.loadings)
        gene_mapping <- setNames(new_gene_names_unique, gene_ids)
        new_loading_names <- gene_mapping[old_loading_names]
        new_loading_names <- new_loading_names[!is.na(new_loading_names)]
        
        if (length(new_loading_names) > 0) {
          # Keep only the loadings that have valid mappings
          valid_indices <- old_loading_names %in% names(gene_mapping)
          reduction@feature.loadings <- reduction@feature.loadings[valid_indices, , drop = FALSE]
          rownames(reduction@feature.loadings) <- new_loading_names
        }
      }
      
      new_seurat@reductions[[reduction_name]] <- reduction
    }
  }
  
  # Update gene annotations in misc to reflect the name change
  if (!is.null(new_seurat@misc$gene_annotations)) {
    # Assuming gene_annotations has the same order as the original genes
    new_seurat@misc$gene_annotations$name <- new_gene_names_unique
    # You might also want to add the original IDs as a separate column
    new_seurat@misc$gene_annotations$original_id <- gene_ids
  }
  
  return(new_seurat)
}

# Usage:
# # Extract the components as you mentioned
# gene_ids <- rownames(sc_org)
# new_gene_names <- sc_org@misc$gene_annotations$name
# 
# # Rename the Seurat object
# sc_org_renamed <- rename_seurat_genes(sc_org, new_gene_names)
# 
# # Verify the rename worked
# print(paste("Original rownames (first 10):", paste(head(gene_ids, 10), collapse = ", ")))
# print(paste("New rownames (first 10):", paste(head(rownames(sc_org_renamed), 10), collapse = ", ")))
# 
# # Check that dimensions are preserved
# print(paste("Original dimensions:", paste(dim(sc_org), collapse = " x ")))
# print(paste("New dimensions:", paste(dim(sc_org_renamed), collapse = " x ")))

# ALTERNATIVE SIMPLER APPROACH (if you have a basic Seurat object):
# This works if you only need to rename the main RNA assay
simple_rename <- function(seurat_obj, new_names) {
  # Make names unique to handle duplicates
  new_names_unique <- make.unique(new_names)
  
  # Rename all matrices in the RNA assay (most common case)
  seurat_obj@assays$RNA@counts@Dimnames[[1]] <- new_names_unique
  seurat_obj@assays$RNA@data@Dimnames[[1]] <- new_names_unique
  
  # Update scale.data if it exists
  if (nrow(seurat_obj@assays$RNA@scale.data) > 0) {
    old_names <- rownames(seurat_obj@assays$RNA@scale.data)
    gene_mapping <- setNames(new_names_unique, rownames(seurat_obj))
    new_scale_names <- gene_mapping[old_names]
    rownames(seurat_obj@assays$RNA@scale.data) <- new_scale_names
  }
  
  # Update meta.features - this is crucial!
  if (!is.null(seurat_obj@assays$RNA@meta.features) && nrow(seurat_obj@assays$RNA@meta.features) > 0) {
    rownames(seurat_obj@assays$RNA@meta.features) <- new_names_unique
  }
  
  # Update variable features
  if (length(VariableFeatures(seurat_obj)) > 0) {
    old_var_features <- VariableFeatures(seurat_obj)
    gene_mapping <- setNames(new_names_unique, rownames(seurat_obj))
    new_var_features <- gene_mapping[old_var_features]
    VariableFeatures(seurat_obj) <- new_var_features[!is.na(new_var_features)]
  }
  
  return(seurat_obj)
}

# Use simple approach:
# sc_org_renamed_simple <- simple_rename(sc_org, new_gene_names)

# QUICK FIX for current error:
# If you're getting the meta.features error right now, run this:
fix_meta_features <- function(seurat_obj) {
  # Get the current rownames from the data matrix
  current_rownames <- rownames(seurat_obj@assays$RNA@data)
  
  # Update meta.features rownames to match
  if (!is.null(seurat_obj@assays$RNA@meta.features) && nrow(seurat_obj@assays$RNA@meta.features) > 0) {
    rownames(seurat_obj@assays$RNA@meta.features) <- current_rownames
  }
  
  return(seurat_obj)
}

# Quick fix usage:
# sc_org <- fix_meta_features(sc_org)
# Then try: sc_org <- NormalizeData(sc_org)