#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#genes to plot expression
genes <- c("MSLN", "CDH1", "CALB2", "MKI67", "SNAI2", "MYBL2", "COL5A1", "AXL", "ITGB4")

#plot monocle umaps
plot_cells(mesoTumor.traj, genes = genes, label_cell_groups = FALSE, show_trajectory_graph = FALSE) + theme_void() +
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
    plot.title = element_text(hjust = 0.5))

