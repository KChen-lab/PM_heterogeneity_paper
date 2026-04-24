#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#Figure4C: umap of Xenium TMA tumor cells only colored by tumor state.

#tumor colors
cell.colors.tumor <- c(
  "Sarcomatoid tumor" = "#006896", "Epithelioid tumor" = "#e84359", "Sarcomatoid proliferating tumor" = "#5abacc", 
  "Epithelioid proliferating tumor" = "#ff8fb3")

#arrow placement on umap
umap <- Embeddings(meso.TMA.tumor, "umap")
xmin <- min(umap[,1])
ymin <- min(umap[,2])

#plot umap
DimPlot(meso.TMA.tumor, group.by = "subcluster_label", label = FALSE, cols = cell.colors.tumor)  +  theme_void() +
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
    plot.title = element_text(hjust = 0.5)) + ggtitle("Tumor Type")
