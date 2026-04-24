#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#figure 3D: umap of scRNA colored by cell type annotations.

#colors
cell.colors <- c(
  "Tumor cells" = "#5abacc", "Alveolar macrophages" = "#afed7c", "NK cells" = "#7B68EE",
  "Treg" = "#6A0DAD", "CD8+ T cells" = "#d62728",  "Interstitial macrophages" = "#32CD32",
  "CD8+ T cells proliferating" = "#9b59b6", "CD4+ T cells" = "#fd6c1d", "Interstitial macrophages M2" = "#8ada0b",
  "Plasma cells" = "#f1c40f", "Interstitial macrophages proliferating" = "#FFEA00", "Classical monocytes" = "#0B6623",
  "Plasmacytoid DC" = "#8c564b", "CD8+ T cells INF active" = "#8B0000", "B cells" = "#bcbd22")

#arrow placement for umap
umap <- Embeddings(meso.scRNA, "umap")
xmin <- min(umap[,1])
ymin <- min(umap[,2])

#plot umap
DimPlot(meso.scRNA, group.by = "refined_clusters", cols = cell.colors) + theme_void() +
  annotate("segment",
           x = xmin, xend = xmin + 2,
           y = ymin, yend = ymin,
           arrow = arrow(length = unit(0.2,"cm"))) +
  annotate("segment",
           x = xmin, xend = xmin,
           y = ymin, yend = ymin + 2,
           arrow = arrow(length = unit(0.2,"cm"))) +
  annotate("text", x = xmin + 2, y = ymin - 0.5, label = "UMAP1") +
  annotate("text", x = xmin - 0.5, y = ymin + 2, label = "UMAP2", angle = 90) + theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)) + ggtitle("Cell Type")
