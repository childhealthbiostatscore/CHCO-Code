# scRNA
```{r echo = F}
# Liver scRNA data processing
# so_liver_sc <- readRDS("/Volumes/Peds Endo/Petter Bjornstad/scRNA/data_raw/PB_Liver_7SingleCellDatasets.RDS")
# so_liver_sc <- readRDS("/Users/hhampson/Dropbox/Bjornstad data/Liver/PB_Liver_7SingleCellDatasets.RDS")
so_liver_sc <- readRDS(fs::path(dir.dat,"scRNA","data_raw","PB_Liver_7SingleCellDatasets.RDS"))

# Add missing liver meta data to Seurat object
# meta_liver_raw <- read.csv("/Volumes/Peds Endo/Petter Bjornstad/scRNA/data_clean/liver_biopsy_metadata_PN.csv")
meta_liver_raw <- read.csv(fs::path(dir.dat,"scRNA","data_clean","liver_biopsy_metadata_PN.csv"))

meta_liver_sc <-  so_liver_sc@meta.data[,1:9] %>%
  dplyr::mutate(Cryostor_ID = gsub("-Liv", "", orig.ident)) %>%
  left_join(meta_liver_raw)

rownames(meta_liver_sc) <- rownames(so_liver_sc@meta.data)

so_liver_sc <- AddMetaData(so_liver_sc, meta_liver_sc)
so_liver_sc <- NormalizeData(so_liver_sc)
so_liver_sc <- ScaleData(so_liver_sc)
# PCA
so_liver_sc <- RunPCA(so_liver_sc, features = VariableFeatures(object = so_liver_sc))
ElbowPlot(so_liver_sc)
# Cluster cells
so_liver_sc <- FindNeighbors(so_liver_sc)
so_liver_sc <- FindClusters(so_liver_sc)
# Perform UMAP and tSNE
so_liver_sc <- RunUMAP(so_liver_sc, dims = 1:15)
DimPlot(so_liver_sc, reduction = "umap") 
```
