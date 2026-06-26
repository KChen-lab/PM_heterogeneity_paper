#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#colors
cell.cycle.colors <- c("G1" = "#00BA38", "S" = "#00BFC4", "G2M" = "#F8766D" , "Undecided" = "#C77CFF")

#plot umap colored by cell cycle state
plot_cells(mesoTumor.traj, label_groups_by_cluster = FALSE, color_cells_by = "celltype_tumor", label_leaves = FALSE,
           label_branch_points = FALSE, label_cell_groups = FALSE) + scale_color_manual(values = cell.cycle.colors)
