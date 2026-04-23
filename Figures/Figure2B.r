#load libraries
set.seed(1234)
library(Seurat)
library(monocle3)
library(tidyverse)

#figure 2B: umap showing MSLN expression in the tumor cell only
#whole slide surgical resection object.

#arrow placement for umap
umap <- reducedDims(meso.WS.Tumor.traj)$UMAP
xmin <- min(umap[,1])
ymin <- min(umap[,2])

#genes to plot
epi.gene <- c("MSLN")

#plot umap with MSLN expression
plot_cells(meso.WS.Tumor.traj, genes = epi.genes, label_cell_groups = FALSE, show_trajectory_graph = FALSE) + theme_void() +
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
    plot.title = element_text(hjust = 0.5)) + ggtitle("MSLN Expression")
