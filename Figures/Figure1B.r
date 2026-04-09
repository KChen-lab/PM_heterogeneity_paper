#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#colors
tumor.colors <- c(
  "Sarcomatoid tumor" = "#006896", "Epithelioid tumor" = "#e84359", "Sarcomatoid proliferating tumor" = "#5abacc", 
  "Epithelioid proliferating tumor" = "#ff8fb3")

#plot UMAP
p <- DimPlot(meso.WS.Tumor, group.by = "refined_tumor_cluster", cols = tumor.colors, label = FALSE, label.size = 4, raster = FALSE) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 5))) +  theme_void() +
  theme(legend.title = element_blank(),  plot.title = element_text(hjust = 0.5)) + ggtitle("Tumor Cell UMAP")

p