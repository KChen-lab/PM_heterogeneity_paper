#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

##panel F: Module scores for the epithelioid, sarcomatoid, and transition regions
#plotted on the whole slide surgical resection visium data.

#create epithelioid, sarcomatoid, and transition signatures
epi.sig <- c("MSLN", "CDH1", "CDH3", "WT1", "CALB2")
sarc.sig <- c("AXL", "COL5A2", "COL5A1", "GAS6", "ITGB4", "LTBP2", "CALB2")
trans.sig <- c("MKI67", "TOP2A", "AURKA", "SNAI2", "TWIST1", "TWIST2", "MYBL2")

#compute module score
vis.WS <- AddModuleScore(vis.WS, features=list(epi.sig), name = "Epithelioid_signature", ctrl = 50)
vis.WS <- AddModuleScore(vis.WS, features=list(sarc.sig), name = "Sarcomatoid_signature", ctrl = 50)
vis.WS<- AddModuleScore(vis.WS, features=list(trans.sig), name = "Transition_signature", ctrl = 50)

#plot module score on visium image
SpatialFeaturePlot(vis.WS, features = c("Epithelioid_signature1", "Sarcomatoid_signature1", "Transition_signature1"), images = "MESO_172T")
