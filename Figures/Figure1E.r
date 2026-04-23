#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

##panel E: Dot plot showing marker gene expression for the four tumor states
#using the tumor cell only whole slide surgical resection object.

#order of cell types on dotplot
cell.type.order <- c("Sarcomatoid tumor", "Sarcomatoid proliferating tumor",
                     "Epithelioid proliferating tumor", "Epithelioid tumor")

#apply order
meso.WS.Tumor$celltype_tumor <- factor(Idents(meso.WS.Tumor), levels = cell.type.order)
Idents(meso.WS.Tumor) <- meso.WS.Tumor$celltype_tumor

#plot dotplot
DotPlot(meso.WS.Tumor, features = c("CDH1", "CDH3", "MSLN", "CALB2", "WT1",
                               "TOP2A", "MKI67", "AURKA", "SNAI2", "TWIST1", "TWIST2", "MYBL2",
                               "THBS2", "COL5A2", "LTBP2", "COL5A1", "AXL", "GAS6", "ITGB4"), 
        group.by="celltype") + RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "#1E90FF") +
  ggplot2::theme( axis.title.x = ggplot2::element_text(size = 8), axis.title.y = ggplot2::element_text(size = 8),
    axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, size = 8), axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::scale_y_discrete(labels = function(x) str_wrap(x, width = 16)) + 
  theme(legend.text = element_text(size = 8)) + theme(legend.title = element_text(size = 8))
