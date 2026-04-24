#load libraries
set.seed(1234)
library(SeuratObject)
library(Seurat)
library(DoubletFinder)
library(tidyverse)

#creat seurat object for the scRNA data from the whole slide surgical resection data (WS).
scRNA.WS <- CreateSeuratObject(counts = "/path/to/counts", min.cells = 10, min.features = 100,
                              assay = "RNA", project = "Whole_slide")
							  
#remove cells with fewer than 200 transcripts, greater that 6000 transcripts, and less than 0.1% mito DNA							  
scRNA.WS <- subset(x = scRNA.WS, subset = nFeature_RNA > 200  & nFeature_RNA < 6000 & percent.mito < 0.1 )

#dimensionality reduction and clustering
scRNA.WS <- scRNA.WS %>% 
  NormalizeData() %<%
  FindVariableFeatures() %>%
  ScaleData() %>%
 
 scRNA.WS <- scRNA.WS %>%
   RunPCA(npcs = 30) %>%
   FindNeighbors(dims = 1:30) %>%
   FindClusters(resolution = 1) %>%
   RunUMAP(dims = 1:30)
   
#find doublets 
#pK Identification (no ground-truth)
sweep.res.list <- paramSweep(scRNA.WS, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

#Homotypic Doublet Proportion Estimate
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075*nrow(scRNA.WS@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#Run DoubletFinder 
scRNA.WS <- doubletFinder(scRNA.WS, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)

#find cell markers per cluster
markers <- FindAllMarkers(scRNA.WS, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

#save object
saveRDS(scRNA.WS, file = "MESO_whole_slide_scRNA.rds")
