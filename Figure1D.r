#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

##panel D
#create tumor subcluster
meso.WS.Tumor <- subset(meso.WS, idents = c("Epithelioid tumor", "Epithelioid proliferating tumor",
                                     "Sarcomatoid proliferating tumor", "Sarcomatoid tumor"))

#crop for transition region
crop <- Crop(meso.WS.Tumor[["fov"]], y = c(8500, 9100), x = c(5200, 6000), coords = "plot")

#apply cropped coords
meso.WS.Tumor[["zoom"]] <- crop
DefaultBoundary(meso.WS.Tumor[["zoom"]]) <- "segmentation"

#cell colors
tumor.colors <- c(
  "Sarcomatoid tumor" = "#006896", "Epithelioid tumor" = "#e84359", "Sarcomatoid proliferating tumor" = "#5abacc", 
  "Epithelioid proliferating tumor" = "#ff8fb3")
  
#plot xenium zoomed image
ImageDimPlot(meso.WS.Tumor, fov = "zoom", group.by = "refined_tumor_cluster", cols = cell.colors, border.size = 0, axes = TRUE)  + 
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
