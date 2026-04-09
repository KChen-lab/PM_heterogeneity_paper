#load libraries
set.seed(1234)
library(cellchat2)
library(ggplot2)
library(tidyverse)

#colors
my_cols <- c(
  "Tumor cells" = "#ADD8E6", "NK cells" = "#7B68EE", "Alveolar macrophages" = "#afed7c",
  "Treg" = "#6A0DAD", "CD8+ T cells" = "#d62728",  "Interstitial mph prevascular" = "#32CD32",
  "Proliferating T cells" = "#9b59b6", "CD4+ T cells" = "#fd6c1d", "Smooth muscle cells" = "#8B0000",
  "CAFs" = "#F39C12", "Plasma cells" = "#f1c40f", "Alveolar mph proliferating" = "#FFEA00",
  "Mast cells" = "#CC5500", "Classical monocytes" = "#0B6623", "EC" = "#C71585",
  "Plasmacytoid DC" = "#8c564b", "DC2" = "#00FF7F",
  "B cells" = "#bcbd22", "Lymphatic EC" = "pink", "Proliferating CAFs" = "#fffb05",
  "Pericytes" = "#CD5C5C")
idents <- levels(cellchat@idents$Biphasic)
my_cols_use <- my_cols[idents]

#plot network graph
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", color.use = my_cols_use,
                          title.name = "Differential Interaction Strength between Epithelioid and Biphasic Samples")
