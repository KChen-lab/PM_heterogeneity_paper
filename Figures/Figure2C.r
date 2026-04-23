#load libraries
set.seed(1234)
library(Seurat)
library(monocle3)
library(tidyverse)

#figure 2C: umap with learned trajectory plotted and colored by pseudotime
#for the tumor cell only whole slide surgical resection xenium object.

#arrow placement for umap
umap <- reducedDims(meso.WS.Tumor.traj)$UMAP
xmin <- min(umap[,1])
ymin <- min(umap[,2])

#plot umap colored by pseudotime
plot_cells(meso.WS.Tumor.traj, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,
           label_branch_points = FALSE) + theme_void() +
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
    plot.title = element_text(hjust = 0.5)) + ggtitle("Pseudotime Trajectory")
