##NEBULA Function 



check_and_load <- function(packages) {
  for(i in c(1:length(packages))){
  if (!requireNamespace(packages[i], quietly = TRUE)) {
    install.packages(packages)  # Install the package if not available
  }
  library(packages[i], character.only = TRUE)  # Load the package
  }
}


libraries <- c('scran', 'future', 'future.apply', 'tidyverse', 
               'colorspace', 'patchwork', 'ggdendro', 'cowplot',
               'ggpubr', 'rstatix', 'arsenal', 'Biobase', 'msigdbr',
               'kableExtra', 'knitr', 'REDCapR', 'data.table', 'emmeans', 
               'NMF', 'pheatmap', 'UpSetR', 'enrichR', 'WriteXLS', 'SAVER',
               'readxl', 'limma', 'edgeR', 'BiocGenerics', 'GSEABase', 'slingshot',
               'SingleCellExperiment', 'MAST', 'muscat', 'scater', 'Seurat', 'jsonlite',
               'dplyr', 'glmmTMB', 'reshape2', 'broom.mixed', 'nebula')




#data is the Seurat object. output is a file path to write to, method is highly variable genes, celltypes can be all or a vector of types. 
#Can choose cell labelling (harmony or rpca)
run_nebula <- function(data, output, offset = 'pooled', method = 'HVG', celltypes = 'All', 
                       cell_labels = 'harmony', hvgcount = 2000){
  check_and_load(libraries)
  if(offset == 'pooled'){
    
    counts_layer <- round(GetAssayData(data, layer = 'counts'))
    library_size <- Matrix::colSums(round(GetAssayData(data, layer = 'counts')))
    data$library_size <- library_size
    
    cat('Calculating pooled offsets. This may take a bit.')
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_layer))
    sce <- computeSumFactors(sce)
    # View size factors
    sizeFactors(sce)
    ## Calculate offset → (size factors)
    data$pooled_offset <- (sizeFactors(sce))
    cat('Done calculating offset. Moving to next analysis')
    
    
  if(celltypes == 'All'){
    
    if(method == 'HVG'){
      cat('Identifying highly variable genes')
    full_analysis <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
    hvgs <- VariableFeatures(full_analysis)
    data<- subset(data, features = hvgs)
    full_counts <- round(GetAssayData(data, layer = "counts")) 
    }
    
    meta_gene <- data@meta.data
    pred_gene <- model.matrix(~group*sex, data = meta_gene)
    data_g_gene <- list(count = data, id = meta_gene$record_id, pred = pred_gene)
    result_allcells <- nebula(count = full_counts, id = data$record_id, pred = data_g_gene$pred, 
                              offset = data$pooled_offset,
                              ncore = 1, output_re = T, covariance = T,
                              reml = T, model = "NBLMM")
    
    
    
    write.table(result_allcells, output, row.names=F,
                quote=F, sep='\t')
    
    return(as.data.frame(result_allcells))
    
    
    
    }else if(celltypes != 'All'){
      
      if(cell_labels == 'harmony'){
        cat('Filtering to specified cell types in Harmony')
      data$celltype_subset <- ifelse(data$celltype_rpca %in% celltypes,
                                      "keep", as.character(data$celltype_rpca))
      data <- subset(data, celltype_ec == "keep")
      }else if(cell_labels == 'rpca'){
        cat('Filtering to specified cell types in RPCA')
        data$celltype_subset <- ifelse(data$celltype_rpca %in% celltypes,
                                       "keep", as.character(data$celltype_rpca))
        data <- subset(data, celltype_ec == "keep")
      }
      
      
      if(method == 'HVG'){
        cat('Identifying highly variable genes')
        full_analysis <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
        hvgs <- VariableFeatures(full_analysis)
        data<- subset(data, features = hvgs)
        full_counts <- round(GetAssayData(data, layer = "counts")) 
      }
      
      meta_gene <- data@meta.data
      pred_gene <- model.matrix(~group*sex, data = meta_gene)
      data_g_gene <- list(count = data, id = meta_gene$record_id, pred = pred_gene)
      result_allcells <- nebula(count = full_counts, id = data$record_id, pred = data_g_gene$pred, 
                                offset = data$pooled_offset,
                                ncore = 1, output_re = T, covariance = T,
                                reml = T, model = "NBLMM")
      
      
      
      write.table(result_allcells, output, row.names=F,
                  quote=F, sep='\t')
      
      return(as.data.frame(result_allcells))      
      
    
    }
    
  }else if(offset == 'library size'){
    cat('Calculating offset by library size instead.')
    counts_layer <- round(GetAssayData(data, layer = 'counts'))
    library_size <- Matrix::colSums(round(GetAssayData(data, layer = 'counts')))
    data$library_size <- library_size
    
    
    if(celltypes == 'All'){
      
      if(method == 'HVG'){
        cat('Identifying highly variable genes')
        full_analysis <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
        hvgs <- VariableFeatures(full_analysis)
        data<- subset(data, features = hvgs)
        full_counts <- round(GetAssayData(data, layer = "counts")) 
      }
      
      meta_gene <- data@meta.data
      pred_gene <- model.matrix(~group*sex, data = meta_gene)
      data_g_gene <- list(count = data, id = meta_gene$record_id, pred = pred_gene)
      result_allcells <- nebula(count = full_counts, id = data$record_id, pred = data_g_gene$pred, 
                                offset = log10(data$library_size),
                                ncore = 1, output_re = T, covariance = T,
                                reml = T, model = "NBLMM")
      
      
      
      write.table(result_allcells, output, row.names=F,
                  quote=F, sep='\t')
      
      return(as.data.frame(result_allcells))
      
      
      
    }else if(celltypes != 'All'){
      
      if(cell_labels == 'harmony'){
        cat('Filtering to specified cell types in Harmony')
        data$celltype_subset <- ifelse(data$celltype_rpca %in% celltypes,
                                       "keep", as.character(data$celltype_rpca))
        data <- subset(data, celltype_ec == "keep")
      }else if(cell_labels == 'rpca'){
        cat('Filtering to specified cell types in RPCA')
        data$celltype_subset <- ifelse(data$celltype_rpca %in% celltypes,
                                       "keep", as.character(data$celltype_rpca))
        data <- subset(data, celltype_ec == "keep")
      }
      
      
      if(method == 'HVG'){
        cat('Identifying highly variable genes')
        full_analysis <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
        hvgs <- VariableFeatures(full_analysis)
        data<- subset(data, features = hvgs)
        full_counts <- round(GetAssayData(data, layer = "counts")) 
      }
      
      meta_gene <- data@meta.data
      pred_gene <- model.matrix(~group*sex, data = meta_gene)
      data_g_gene <- list(count = data, id = meta_gene$record_id, pred = pred_gene)
      result_allcells <- nebula(count = full_counts, id = data$record_id, pred = data_g_gene$pred, 
                                offset = log10(data$library_size),
                                ncore = 1, output_re = T, covariance = T,
                                reml = T, model = "NBLMM")
      
      
      
      write.table(result_allcells, output, row.names=F,
                  quote=F, sep='\t')
      
      return(as.data.frame(result_allcells))      
      
      
    }
    
    
    
    
    
  }
  
  
}




















