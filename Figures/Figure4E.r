#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#p1.3
min.x = 5830.020
max.x = 7101.974
min.y = 5812.758
max.y = 7103.752

#crop image 
crop <- Crop(meso[["fov"]], x = c(min.x, max.x), y = c(min.y, max.y), coords = "tissue")

#add cropped image to new assay plot
p1.3 <- meso
p1.3[["fov"]] <- crop
DefaultBoundary(p1.3[["fov"]]) <- "segmentation"

#subset cells within the crop coordinates
coords <- GetTissueCoordinates(p1.3)
in.bounds <- coords$x >= min.x & coords$x <= max.x &
  coords$y >= min.y & coords$y <= max.y
coords.in.bounds <- coords[in.bounds, ]
cells.in.crop <- coords.in.bounds$cell

#subset seurat object by cells in crop
p1.3 <- subset(p1.3, cells = cells.in.crop)

#cell colors
cell.colors <- c(
  "Sarcomatoid tumor" = "#006896", "Epithelioid tumor" = "#e84359", "Sarcomatoid proliferating tumor" = "#5abacc", 
  "Epithelioid proliferating tumor" = "#ff8fb3", "NK cells" = "#7B68EE", "Alveolar macrophages" = "#afed7c",
  "Treg" = "#6A0DAD", "CD8+ T cells" = "#d62728",  "Interstitial mph prevascular" = "#32CD32",
  "Proliferating CD8+ T cells" = "#9b59b6", "CD4+ T cells" = "#fd6c1d", "Smooth muscle cells" = "#8B0000",
  "CAFs" = "#F39C12", "Plasma cells" = "#f1c40f", "Alveolar mph proliferating" = "#FFEA00",
  "Mast cells" = "#CC5500", "Classical monocytes" = "#0B6623", "EC" = "#C71585", 
  "Plasmacytoid DC" = "#8c564b", "DC2" = "#00FF7F", "Proliferating CAFs" = "#fffb05",
  "B cells" = "#bcbd22", "Lymphatic EC" = "pink", "Pericytes" = "#CD5C5C",
  "Ambiguous" = "#E6E6FE")

#check cropped coordinates
ImageDimPlot(p1.3, fov = "fov", axes = TRUE, border.size = 0, cols = cell.colors,
             group.by = "refined_tumor_cluster", coord.fixed = TRUE) + theme(panel.grid = element_blank())

