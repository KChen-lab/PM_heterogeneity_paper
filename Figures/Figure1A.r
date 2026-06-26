#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#panel A: umap of all cell types annotated in the whole slide surgical
#resection Visium spot level data. vis.WS is the visium object generated in
#the whole_slide_visium_object.r file in the Seurat_objects_folder.

#cell type colors
cell.colors <- c(
  "Sarcomatoid Tumor" = "#006896", "Epithelioid Tumor" = "#e84359", "Tumor Transition" = "#5abacc",
  "Interstitial Macrophages" = "#32CD32", "Smooth Muscle cells" = "#8B0000",
  "Fibroblasts" = "#F39C12", "Stromal cells" = "#f1c40f", "B/Plasma cells" = "#bcbd22",
  "Sarcomatoid-Immune Tumor" = "#afed7c")

#plot UMAP
DimPlot(vis.WS, group.by = "refined_tumor_cluster", cols = cell.colors, label = FALSE, label.size = 4, raster = FALSE) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 5))) +  theme_void() + 
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) + ggtitle("Cell Type UMAP")
  
  





