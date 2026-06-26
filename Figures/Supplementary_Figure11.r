#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#colors
patient.colors <- c(
  "1" = "#ff3155", "2" = "#ffaf42", "3" = "#ffed5e", "4" = "#49f770", "5" = "#2daefd", "6" = "#9d4dff")

#arrow placement on umap
umap <- Embeddings(meso, "umap")
xmin <- min(umap[,1])
ymin <- min(umap[,2])

#umap for patient ID
DimPlot(meso, group.by = "patient_number", cols = patient.colors, label = FALSE, label.size = 4, raster = FALSE) +
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
    plot.title = element_text(hjust = 0.5)) + ggtitle("Cell Type Annotation")
