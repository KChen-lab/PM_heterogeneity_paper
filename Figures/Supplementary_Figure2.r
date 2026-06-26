#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#plot dotplot
DotPlot(vis.WS, features = c("TGFB1", "TGFB2", "TGFB3", "MSLN", "WT1", "CALB2", "CLDN15", "CDH1", "CDH3",
                           "TOP2A", "MKI67", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "MYBL2", "AXL", "GAS6", "ITGB4", "GZMB", "CXCL10", "CD4", "LAG3", "MARCO", "AGER", "VWF",
                           "LYZ", "CD163", "CD68", "LYVE1", "SLIT3", "FBN1", "ACTA2", "FBLN2", "MS4A1", "MZB1",
                           "CD79A", "BANK1", "JCHAIN", "PDGFRB", "MFAP5", "MFAP4", "FAP", "MMP3", "MYH11", "CNN1", "TAGLN",
                           "NOTCH3"),
        group.by = "refined_tumor_cluster") +
  scale_color_gradient(low = "lightgrey", high = "#1E90FF") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) +
  ggplot2::scale_y_discrete(labels = function(x) str_wrap(x, width = 16))

