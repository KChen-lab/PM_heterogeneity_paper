#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#dot plot for tumor cells
DotPlot(meso.scRNA.tumor, features = c("MSLN", "CALB2","WT1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", 
                                 "AXL", "GAS6", "COL5A1", "ITGB4"), group.by= "phenotype") + RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "#1E90FF") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, size = 10, hjust = 1))
