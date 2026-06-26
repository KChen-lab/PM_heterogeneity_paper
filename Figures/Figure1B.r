#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#figure 1B:  Spaitial image of cell type predictions from the whole slide
#surgical resection visium spot level data. vis.WS is the 

#colors
cell.colors <- c(
  "Sarcomatoid Tumor" = "#006896", "Epithelioid Tumor" = "#e84359", "Tumor Transition" = "#5abacc",
  "Interstitial Macrophages" = "#32CD32", "Smooth Muscle cells" = "#8B0000",
  "Fibroblasts" = "#F39C12", "Stromal cells" = "#f1c40f", "B/Plasma cells" = "#bcbd22",
  "Sarcomatoid-Immune Tumor" = "#afed7c")

#plot spatial image
SpatialDimPlot(vis.WS, group.by = "refined_tumor_cluster", cols = cell.colors)
