#load libraries
set.seed(1234)
library(Seurat)
library(monocle3)
library(tidyverse)

#figure 2A: umap of tumor cells states with the learned trajectory graph plotted
#on top, created using monocle3 on the tumor cells only whole slide surgical resection object.

#arrow placement for umap
umap <- reducedDims(meso.WS.Tumor.traj)$UMAP
xmin <- min(umap[,1])
ymin <- min(umap[,2])

#colors
tumor.colors <- c(
  "Sarcomatoid tumor" = "#006896", "Epithelioid tumor" = "#e84359", "Sarcomatoid proliferating tumor" = "#5abacc", 
  "Epithelioid proliferating tumor" = "#ff8fb3")

#plot learned trajectory umap
plot_cells(meso.WS.Tumor.traj, color_cells_by = "celltype_tumor", label_groups_by_cluster = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE) + theme_void() +
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
    plot.title = element_text(hjust = 0.5)) + ggtitle("Learned Trajectory Graph") + scale_color_manual(values = tumor.colors)
