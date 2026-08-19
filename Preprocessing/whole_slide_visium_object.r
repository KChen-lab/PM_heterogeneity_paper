#load libraries
set.seed(1234)
library(SeuratObject)
library(Seurat)
library(tidyverse)

#create seurat object from visium whole slide (WS)
vis.WS <- LoadData ("/path to data/",)

#subset spatial features to include >200 features per spot
vis.WS <- vis.WS[, vis.WS$n_Feature_Spatial > 200]

#normalization, dimensionalty reduction, and clustering
vis.WS <- vis.WS %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

vis.WS <- vis.WS %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resoltuion = 0.8) %>%
  RunUMAP(dims = 1:30)

#find cell markers per cluster
Idents(vis.WS) <- "seurat_clusters"
markers <- FindAllMarkers(vis.WS, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

#save object
saveRDS(vis.WS, file = "MESO_Whole_Slide_visium.rds")
