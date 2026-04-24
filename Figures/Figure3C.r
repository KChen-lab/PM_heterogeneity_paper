#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#Figure3C: umap of scRNA seq colored by numbat clones.

#clone colors
clone.colors <- c(
  "1" = "#55d6be", "2" = "#a29bfe", "3" = "#4e6cff")

#arrow placement for umap
umap <- Embeddings(meso.scRNA, "umap")
xmin <- min(umap[,1])
ymin <- min(umap[,2])

#plot umap colored by clone
DimPlot(meso.scRNA, group.by = "clones", cols = clone.colors, cells = colnames(meso.scRNA)[!is.na(meso.scRNA$clones)]) + 
  theme_void() +
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
    plot.title = element_text(hjust = 0.5)) + ggtitle("Clone") + theme(legend.position = "none")
