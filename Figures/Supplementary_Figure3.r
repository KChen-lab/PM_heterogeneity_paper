#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#cell type colors
cell.colors <- c(
  "Tumor cells" = "#5abacc", "NK cells" = "#7B68EE", "Alveolar macrophages" = "#afed7c",
  "Treg" = "#6A0DAD", "CD8+ T cells" = "#d62728",  "Interstitial mph prevascular" = "#32CD32",
  "Proliferating CD8+ T cells" = "#9b59b6", "CD4+ T cells" = "#fd6c1d", "Smooth muscle cells" = "#8B0000",
  "CAFs" = "#F39C12", "Plasma cells" = "#f1c40f", "Alveolar mph proliferating" = "#FFEA00",
  "Mast cells" = "#CC5500", "Classical monocytes" = "#0B6623", "EC" = "#C71585", 
  "Plasmacytoid DC" = "#8c564b", "DC2" = "#00FF7F",
  "B cells" = "#bcbd22", "Lymphatic EC" = "pink", "Pericytes" = "#CD5C5C",
  "Ambiguous" = "#E6E6FE", "Proliferating Treg" = "#1E90FF", "Proliferating NK cells" = "magenta")

#collapse tumor annotations
meso.WS$tumor_collapsed <- meso.WS$refined_tumor_cluster
tumor.levels <- c(
  "Epithelioid tumor",
  "Epithelioid proliferating tumor",
  "Sarcomatoid proliferating tumor",
  "Sarcomatoid tumor")

meso.WS$tumor_collapsed[
  meso.WS$tumor_collapsed %in% tumor.levels
] <- "Tumor cells"

#plot umap
DimPlot(meso.WS, group.by = "tumor_collapsed", cols = cell.colors, label = FALSE, label.size = 4, raster = FALSE) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 5)))
