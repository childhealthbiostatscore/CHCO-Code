#Single Cell RNA Sequencing Quality Control & Pre-Processing Function ----
qc_function <- function(so,var_pct,gene_pct) {
# qc_function <- function(so,cut_low,cut_high,mt_pct,var_pct,gene_pct) {
  
#Find the percent of RNA per cell that is mitochondrial RNA
# so[['percent_mt']] <- PercentageFeatureSet(so, pattern = '^MT-')

#Filter out cells with fewer than "cut_low" genes, cells with "cut_high" detected genes, and cells where < mt_pct% of the expressed genes are mitochondrial genes
# so <- subset(so, subset = nFeature_RNA > cut_low & nFeature_RNA < cut_high & percent_mt<mt_pct)  

#Filter out rare genes expressed in less than "gene_pct" of cells
expr_matrix <- as.matrix(GetAssayData(so, layer = "counts"))
# Calculate the proportion of cells expressing each gene
num_cells_per_gene <- rowSums(expr_matrix > 0)  # Count nonzero cells per gene
total_cells <- ncol(expr_matrix)  # Total number of cells
gene_proportion <- num_cells_per_gene / total_cells  # Fraction of cells expressing each gene
remove(expr_matrix)
# Keep genes expressed in at least "gene_pct" of cells
genes_to_keep <- names(gene_proportion[gene_proportion >= gene_pct])
so <- subset(so, features = genes_to_keep)

#Remove Mitochondrial Genes
mito_genes <- grep("^MT-", rownames(so), value = TRUE)
so <- subset(so, features = setdiff(rownames(so@assays$RNA@layers$counts), mito_genes))

# Normalize the adjusted counts (considered adjusted counts after SoupX processing)
so <- NormalizeData(so)

#Scale the adjusted counts (subtract the mean gene expression of the gene across all cells from the gene adjusted count value, then divide by SD of the gene across all cells) - each gene should have mean of 0, SD of 1
so <- ScaleData(so, features = VariableFeatures(so))

#Perform PCA to reduce the complexity of the data & capture the most important sources of variation
so <- RunPCA(so, features = VariableFeatures(so))
#Get the standard deviation for each Principle component
pca_results <- so[["pca"]]@stdev  # Standard deviation for each PC
#Calculate the proportion of variation explained by each PC
var_explained <- pca_results^2 / sum(pca_results^2)  # Proportion of variance explained
# Choose the number of PCs for "var_pct"% cummulative variance explained
cumulative_var_explained <- cumsum(var_explained)
dims <- which(cumulative_var_explained >= var_pct)[1]  # Select the number of PCs explaining 90% of variance
elbow_plot <- ElbowPlot(so) + ggtitle(label=paste0(dims, " Principle Components explain, ",var_pct*100,"% of Cummulative Variance"))
#Compute the nearest neighbors of each cell based on PCA
so <- FindNeighbors(so, dims = 1:dims)
#Clusters cells that are similar to each other based on their gene expression profiles using nearest neighbors information 
so <- FindClusters(so)
#Perform Uniform Manifold Approximation and Projection (UMAP) - So that you can visualize high dim data on a low dim space (UMAP) using cluster or pca info 
so <- RunUMAP(so, dims = 1:dims, reduction.key = "UMAP_")

# return(so)
return(list(so = so, elbow_plot = elbow_plot))
}

#Single Nuc RNA Sequencing Quality Control & Pre-Processing Function ----
qc_function_sn <- function(so,var_pct) {
  # qc_function <- function(so,cut_low,cut_high,mt_pct,var_pct,gene_pct) {
  
  #Find the percent of RNA per cell that is mitochondrial RNA
  # so[['percent_mt']] <- PercentageFeatureSet(so, pattern = '^MT-')
  
  #Filter out cells with fewer than "cut_low" genes, cells with "cut_high" detected genes, and cells where < mt_pct% of the expressed genes are mitochondrial genes
  # so <- subset(so, subset = nFeature_RNA > cut_low & nFeature_RNA < cut_high & percent_mt<mt_pct)  
  
  # #Filter out rare genes expressed in less than "gene_pct" of cells
  # expr_matrix <- as.matrix(GetAssayData(so, layer = "data"))
  # # Calculate the proportion of cells expressing each gene
  # num_cells_per_gene <- rowSums(expr_matrix > 0)  # Count nonzero cells per gene
  # total_cells <- ncol(expr_matrix)  # Total number of cells
  # gene_proportion <- num_cells_per_gene / total_cells  # Fraction of cells expressing each gene
  # remove(expr_matrix)
  # # Keep genes expressed in at least "gene_pct" of cells
  # genes_to_keep <- names(gene_proportion[gene_proportion >= gene_pct])
  # so <- subset(so, features = genes_to_keep)
  
  #Remove Mitochondrial Genes
  mito_genes <- grep("^MT-", rownames(so), value = TRUE)
  so <- subset(so, features = setdiff(rownames(so@assays$RNA@data), mito_genes))
  
  # Normalize the adjusted counts (considered adjusted counts after SoupX processing)
  so <- NormalizeData(so)
  
  #Scale the adjusted counts (subtract the mean gene expression of the gene across all cells from the gene adjusted count value, then divide by SD of the gene across all cells) - each gene should have mean of 0, SD of 1
  so <- ScaleData(so, features = VariableFeatures(so))
  
  #Perform PCA to reduce the complexity of the data & capture the most important sources of variation
  so <- RunPCA(so, features = VariableFeatures(so))
  #Get the standard deviation for each Principle component
  pca_results <- so[["pca"]]@stdev  # Standard deviation for each PC
  #Calculate the proportion of variation explained by each PC
  var_explained <- pca_results^2 / sum(pca_results^2)  # Proportion of variance explained
  # Choose the number of PCs for "var_pct"% cummulative variance explained
  cumulative_var_explained <- cumsum(var_explained)
  dims <- which(cumulative_var_explained >= var_pct)[1]  # Select the number of PCs explaining 90% of variance
  elbow_plot <- ElbowPlot(so) + ggtitle(label=paste0(dims, " Principle Components explain, ",var_pct*100,"% of Cummulative Variance"))
  #Compute the nearest neighbors of each cell based on PCA
  so <- FindNeighbors(so, dims = 1:dims)
  #Clusters cells that are similar to each other based on their gene expression profiles using nearest neighbors information 
  so <- FindClusters(so)
  #Perform Uniform Manifold Approximation and Projection (UMAP) - So that you can visualize high dim data on a low dim space (UMAP) using cluster or pca info 
  so <- RunUMAP(so, dims = 1:dims, reduction.key = "UMAP_")
  
  # return(so)
  return(list(so = so, elbow_plot = elbow_plot))
}


#UMAP Visualization Function ----
umap_vis <- function(so,cell_type,group_variable) {
dimplot1 <- DimPlot(so, reduction = "umap",group.by = cell_type,label=T,raster=F) + 
  ggtitle(paste0("UMAP by ",str_to_title(cell_type)))

dimplot2 <- DimPlot(so, reduction = "umap",group.by = group_variable,label=F,raster=F) + 
  ggtitle(paste0("UMAP by ",str_to_title(group_variable)))

#Print out all visualizations
print(dimplot1/dimplot2)
}

#Gene Expression Extraction Function ----
#Extract gene expression matrix to perform linear mixed effect models, etc. 
##a. All cell types ----
ext_gene_fxn <- function(so) {
  for (celltype in so_filtered$KPMP_celltype) {
    so_celltype <- subset(so_filtered,KPMP_celltype %in% celltype)
    
    # Extract the gene expression data for all genes
    # gene_expression <- as.data.frame(so_celltype@assays$RNA@layers$counts)
    gene_expression <- as.data.frame(GetAssayData(so_filtered, layer = "data"))
    
    #Assign gene names as rownames, cell names as colnames
    # rownames(gene_expression) <- rownames(so_celltype) #Gene Names
    # colnames(gene_expression) <- colnames(so_celltype) #Cell Names
    #Transpose gene expression dataset to merge with metadata
    gene_expression <- t(gene_expression)
    gene_expression <- data.frame(gene_expression) #Make a dataframe again after transposing
    #Set gene list
    gene_list <- colnames(gene_expression)
    gene_expression$cellname <- rownames(gene_expression) 
    rownames(gene_expression) <- NULL
    
    # Extract the metadata
    metadata <- so_celltype@meta.data
    # metadata <- metadata %>%
    #   mutate(across(everything(),~ifelse(.==".",NA,.)))
    metadata$cellname <- rownames(metadata)
    rownames(metadata) <- NULL
    
    # Combine the gene expression data and metadata
    data <- tidylog::full_join(metadata,gene_expression,by="cellname")
    rm(metadata,gene_expression)
    
    N <- length(unique(data$kit_id))
  }
}
##b. By Specific Cell Type ----
#Must change the cell type variable in the function if using other datasets where celltype variable isnt called KPMP_celltype
ext_gene_fxn2 <- function(so,celltype) {
# for (celltype in so_filtered$KPMP_celltype) {
  # so_celltype <- subset(so,KPMP_celltype %in% celltype)
  so_celltype <- subset(so,celltype2 %in% celltype)
  
  
  # Extract the gene expression data for all genes
  gene_expression <- as.data.frame(GetAssayData(so_celltype, layer = "data"))
  
  #Assign gene names as rownames, cell names as colnames
  # rownames(gene_expression) <- rownames(so_celltype) #Gene Names
  # colnames(gene_expression) <- colnames(so_celltype) #Cell Names
  #Transpose gene expression dataset to merge with metadata
  gene_expression <- t(gene_expression)
  gene_expression <- data.frame(gene_expression) #Make a dataframe again after transposing
  #Set gene list
  gene_list <- colnames(gene_expression)
  gene_expression$cellname <- rownames(gene_expression) 
  rownames(gene_expression) <- NULL
  
  # Extract the metadata
  metadata <- so_celltype@meta.data
  # metadata <- metadata %>%
  #   mutate(across(everything(),~ifelse(.==".",NA,.)))
  metadata$cellname <- rownames(metadata)
  rownames(metadata) <- NULL
  
  # Combine the gene expression data and metadata
  data <- tidylog::full_join(metadata,gene_expression,by="cellname")
  rm(metadata,gene_expression)
  
  # N <- length(unique(data$kit_id))
  return(list(data=data, gene_list=gene_list))
# }
}

#Linear Mixed-Effects Model Function ----
# lmm_function <- function()