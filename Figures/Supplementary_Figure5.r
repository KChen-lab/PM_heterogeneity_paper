#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#cell colors
tumor.colors <- c(
  "Sarcomatoid tumor" = "#006896", "Epithelioid tumor" = "#e84359", "Sarcomatoid proliferating tumor" = "#5abacc", 
  "Epithelioid proliferating tumor" = "#ff8fb3")

#plot umap
DimPlot(meso.WS.Tumor, group.by = "refined_tumor_annotation", cols = tumor.colors, label = FALSE, label.size = 4, raster = FALSE) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 5)))


