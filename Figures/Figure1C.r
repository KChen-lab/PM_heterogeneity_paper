#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

##panel C
#cell type colors
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

#plot xenium cell type annotated image
ImageDimPlot(meso.WS, fov = "fov", group.by = "refined_tumor_cluster", cols = cell.colors, border.size = 0) + 
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) + ggtitle("Cell Type Annotation")
