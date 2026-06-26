#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#Supplementary Figure 12: whole slide surgical resection Xenium imaged cropped to tissue core hole area.
#Single core TMA image corresponding to tissue core whole. Colored by cell type annnotation.

#crop for tissue core hole
crop <- Crop(meso.WS[["fov"]], y = c(7000, 9100), x = c(2000, 4000), coords = "plot")

#apply cropped coords
meso.WS[["zoom"]] <- crop
DefaultBoundary(meso.WS[["zoom"]]) <- "segmentation"

#colors
cell.colors <- c(
  "Sarcomatoid tumor" = "#006896", "Epithelioid tumor" = "#e84359", "Sarcomatoid proliferating tumor" = "#5abacc", 
  "Epithelioid proliferating tumor" = "#ff8fb3", "NK cells" = "#7B68EE", "Alveolar macrophages" = "#afed7c",
  "Treg" = "#6A0DAD", "CD8+ T cells" = "#d62728",  "Interstitial mph prevascular" = "#32CD32",
  "Proliferating CD8+ T cells" = "#9b59b6", "CD4+ T cells" = "#fd6c1d", "Smooth muscle cells" = "#8B0000",
  "CAFs" = "#F39C12", "Plasma cells" = "#f1c40f", "Alveolar mph proliferating" = "#FFEA00",
  "Mast cells" = "#CC5500", "Classical monocytes" = "#0B6623", "EC" = "#C71585", 
  "Plasmacytoid DC" = "#8c564b", "DC2" = "#00FF7F",
  "B cells" = "#bcbd22", "Lymphatic EC" = "pink", "Pericytes" = "#CD5C5C",
  "Ambiguous" = "#E6E6FE", "Proliferating Treg" = "#1E90FF", "Proliferating NK cells" = "magenta")

#image plot
ImageDimPlot(meso.WS, fov = "zoom", group.by = "refined_tumor_cluster", cols = cell.colors, border.size = 0.05, na.value = "black", border.color = "black") +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) + theme(legend.position = "none")

#crop for tissue core
#p6.1
min.x = 8855.367
max.x = 10112.088
min.y = 1368.732
max.y = 2667.343

#crop image 
crop <- Crop(meso.TMA[["fov"]], x = c(min.x, max.x), y = c(min.y, max.y), coords = "tissue")

#add cropped image to new assay plot
p6.1 <- meso.TMA
p6.1[["fov"]] <- crop
DefaultBoundary(p6.1[["fov"]]) <- "segmentation"

#subset cells within the crop coordinates
coords <- GetTissueCoordinates(p6.1)
in.bounds <- coords$x >= min.x & coords$x <= max.x &
  coords$y >= min.y & coords$y <= max.y
coords.in.bounds <- coords[in.bounds, ]
cells.in.crop <- coords.in.bounds$cell

#subset seurat object by cells in crop
p6.1 <- subset(p6.1, cells = cells.in.crop)

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
ImageDimPlot(p6.1, fov = "fov", axes = TRUE, border.size = 0, cols = cell.colors,
             group.by = "refined_tumor_cluster", coord.fixed = TRUE) + theme(panel.grid = element_blank())
