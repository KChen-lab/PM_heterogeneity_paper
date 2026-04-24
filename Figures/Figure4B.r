#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#Figure 4B: umap of Xenium TMA data colored by cell type annotation.

#colors
cell.colors <- c(
  "Tumor cells" = "#5abacc", "NK cells" = "#7B68EE", "Alveolar macrophages" = "#afed7c",
  "Treg" = "#6A0DAD", "CD8+ T cells" = "#d62728",  "Interstitial mph prevascular" = "#32CD32",
  "Proliferating T cells" = "#9b59b6", "CD4+ T cells" = "#fd6c1d", "Smooth muscle cells" = "#8B0000",
  "CAFs" = "#F39C12", "Plasma cells" = "#f1c40f", "Alveolar mph proliferating" = "#FFEA00",
  "Mast cells" = "#CC5500", "Classical monocytes" = "#0B6623", "EC" = "#C71585",
  "Plasmacytoid DC" = "#8c564b", "DC2" = "#00FF7F", "Nonclassical monocytes" = "#0000CD",
  "B cells" = "#bcbd22", "Lymphatic EC" = "pink", "Proliferating CAFs" = "#fffb05",
  "Pericytes" = "#CD5C5C", "Ambiguous" = "#E6E6FE")

#arrow placement on UMAP
umap <- Embeddings(meso.TMA, "umap")
xmin <- min(umap[,1])
ymin <- min(umap[,2])

#plot UMAP
DimPlot(meso.TMA, group.by = "tumor_collapsed", cols = cell.colors, label = FALSE, label.size = 4, raster = FALSE) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 5))) +  theme_void() +
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
    plot.title = element_text(hjust = 0.5)) + ggtitle("Cell Type Annotation") + theme(legend.position = "none")
