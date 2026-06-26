#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#plot dot plot
DotPlot(scRNA.WS, features = c("MSLN", "CALB2","WT1", "AXL", "GAS6","COL5A1", "ITGB4",
                                 "CD14", "CD163", "MRC1", "MARCO", "TREM2", "FCGR1A", "ETV5", "CD68", "LILRB2",
                                 "ITGAM", "MS4A1", "CD19", "CD79A", "MZB1", "IRF4",
                                 "CD38","NCAM1", "KLRD1", "NKG7", "CD2", "CD3E", "CD4", "CD8A",
                                 "FOXP3", "CTLA4", "TIGIT","TOP2A", "MKI67"), group.by= "refined_clusters") + RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "#1E90FF") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, size = 10, hjust = 1))
