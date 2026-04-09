#load libraries
set.seed(1234)
library(SeuratObject)
library(Seurat)
library(tidyverse)

#create seurat object for the tissue microarray (TMA)
meso.TMA) <- LoadXenium("/path to data/", fov = "fov", assay = "Xenium")

#remove cells with fewer than 10 transcript counts
meso.TMA <- subset(meso.TMA, subset = nCount_Xenium > 10)

#set default assay and boundary
DefaultAssay(meso.TMA) <- "Xenium"
DefaultBoundary(meso.TMA[["fov"]]) <- "segmentation"

#median transcript counts per cell
counts <- GetAssayData(meso.TMA, assay = "Xenium", layer = "counts")
cell.counts <- Matrix::colSums(cts)
median(cell.counts)

#normalization, dimensionalty reduction, and clustering
meso.TMA <- meso.TMA %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

meso.TMA <- meso.TMA %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resoltuion = 0.8) %>%
  RunUMAP(dims = 1:30)
  
#find cell markers per cluster 
markers <- FindAllMarkers(meso.TMA, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

#save object
saveRDS(meso.TMA, file = "MESO_TMA_Xenium.rds")