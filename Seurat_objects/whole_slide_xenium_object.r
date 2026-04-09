#load libraries
set.seed(1234)
library(SeuratObject)
library(Seurat)
library(tidyverse)

#create seurat object for the whole slide (WS)
meso.WS <- LoadXenium("/path to data/", fov = "fov", assay = "Xenium")

#remove cells with fewer than 10 transcript counts
meso.WS <- subset(meso.WS, subset = nCount_Xenium > 10)

#set default assay and boundary
DefaultAssay(meso.WS) <- "Xenium"
DefaultBoundary(meso.WS[["fov"]]) <- "segmentation"

#median transcript counts per cell
counts <- GetAssayData(meso.WS, assay = "Xenium", layer = "counts")
cell.counts <- Matrix::colSums(cts)
median(cell.counts)

#normalization, dimensionalty reduction, and clustering
meso.WS <- meso.WS %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

meso.WS <- meso.WS %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resoltuion = 0.8) %>%
  RunUMAP(dims = 1:30)
  
#find cell markers per cluster 
markers <- FindAllMarkers(meso.WS, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

#save object
saveRDS(meso.WS, file = "MESO_whole_slide_Xenium.rds")