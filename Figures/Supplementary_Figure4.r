#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#genes to plot 
features.to.plot <- c("CLDN5", "CD34", "ADGRL4", "PDGFRB", "CSPG4", "RGS5", "SLIT3", 
                      "PDGFRA", "POSTN", "PDPN", "FAP", "MMP11", "COMP",
                      "THBS2", "MYBL2", "LTBP2", "AXL", "GAS6", "ITGB4",
                      "COL5A1", "COL5A2", "SNAI2", "TWIST1", "TWIST2", "TOP2A", "MKI67",
                      "MSLN", "CALB2", "WT1", "CDH1", "CDH3",
                      "CD14", "CD163", "ETV5", "CD68", "LILRB2",
                      "ITGAM", "ITGAX", "MS4A1", "CD79A", "MZB1", "IRF4", "CD19", "CD1C", "CLEC10A",
                      "KIT", "MS4A2", "CD38", "KLRD1", "CD2", "CD3E", "CD4", "CD8A", "CD247",
                      "FOXP3", "CTLA4", "TIGIT")

cell.type.order <- c("Ambiguous", "Proliferating NK cells", "NK cells", "Proliferating Treg", "Treg", "CD4+ T cells", 
                     "Proliferating CD8+ T cells", "CD8+ T cells", "Mast cells", "DC2", "Plasmacytoid DC", 
                     "Plasma cells", "B cells", "Classical monocytes", "Alveolar mph proliferating", 
                     "Alveolar macrophages", "Interstitial mph prevascular", "Epithelioid tumor",
                     "Epithelioid proliferating tumor", "Sarcomatoid proliferating tumor", "Sarcomatoid tumor", 
                     "CAFs", "Smooth muscle cells", "Pericytes", "Lymphatic EC", "EC")

#apply order
meso.WS$celltype_tumor <- factor(Idents(meso.WS), levels = cell.type.order)
Idents(meso.WS) <- meso.WS$celltype_tumor

#plot dotplot
DotPlot(meso.WS, features = features.to.plot) +
  scale_color_gradient(low = "lightgrey", high = "#1E90FF") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, size = 8, hjust = 1))
  
