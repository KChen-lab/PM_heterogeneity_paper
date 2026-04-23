#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(viridis)

#Figure 2E: four cropped xenium images colored by pseudotime bins.

#crop for morphology transition zoomed
crop1 <- Crop(meso.WS.Tumor[["fov"]], y = c(8600, 9000), x = c(5300, 5800), coords = "plot")
#apply cropped coords
meso.WS.Tumor[["zoom1"]] <- crop1
DefaultBoundary(meso.WS.Tumor[["zoom1"]]) <- "segmentation"

#crop for morphology transition zoomed epithelioid
crop2 <- Crop(meso.WS.Tumor[["fov"]], y = c(8600, 8800), x = c(5300, 5600), coords = "plot")
#apply cropped coords
meso.WS.Tumor[["zoom2"]] <- crop2
DefaultBoundary(meso.WS.Tumor[["zoom2"]]) <- "segmentation"

#crop for morphology transition zoomed transition 
crop3 <- Crop(meso.WS.Tumor[["fov"]], y = c(8700, 8900), x = c(5400, 5700), coords = "plot")
#apply cropped coords
meso.WS.Tumor[["zoom3"]] <- crop3
DefaultBoundary(meso.WS.Tumor[["zoom3"]]) <- "segmentation"

#crop for morphology transition zoomed sarcomatoid
crop4 <- Crop(meso.WS.Tumor[["fov"]], y = c(8800, 9000), x = c(5500, 5800), coords = "plot")
#apply cropped coords
meso.WS.Tumor[["zoom4"]] <- crop4
DefaultBoundary(meso.WS.Tumor[["zoom4"]]) <- "segmentation"

#plot colors
plasma.grad <- plasma(40)
purple.hue <- plasma.grad[5]
pink.hue <- plasma.grad[18]
yellow.hue <- plasma.grad[33]

##plot crop1
#colors
tumor.colors <- c(
  "Sarcomatoid tumor" = "#006896", "Epithelioid tumor" = "#e84359", "Sarcomatoid proliferating tumor" = "#5abacc", 
  "Epithelioid proliferating tumor" = "#ff8fb3")
  
#plot with correct segmentation boundaries
p1 <- ImageDimPlot(meso.WS.Tumor, fov = "zoom1", group.by = "refined_tumor_cluster", cols = tumor.colors, border.size = 0, na.value = "grey80") +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())

#plot with incorrect segmentation boundaries
p2 <- ImageFeaturePlot(meso.WS.Tumor, fov = "zoom1", features = "pseudotime", cols = plasma.grad, dark.background = TRUE, 
                       border.size = 0.15, boundaries = "segmentations", border.color = "black") + 
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) + 
  theme(legend.position = "none") + labs(title = NULL)

#add p1 (x,y) coordinates to p2
p2[[1]]$data <- p2[[1]]$data %>%
  left_join(p1[[1]]$data %>% select(cell, x_correct = x, y_correct = y),
            by = "cell") %>%
  mutate(x = x_correct, y = y_correct,
         x_correct = NULL, y_correct = NULL)

p2

#plot crop2
ImageDimPlot(meso.WS.Tumor, fov = "zoom2", group.by = "bin_Epithelioid", cols = c("gray80", purple.hue), border.size = 0.1,
  border.color = "white", axes = FALSE, na.value = "grey80") + theme(legend.position = "none") + labs(title = NULL)
  
 #plot crop3
ImageDimPlot(meso.WS.Tumor, fov = "zoom3", group.by = "bin_Transition", cols = c("gray80", pink.hue), border.size = 0.1,
  border.color = "white", axes = FALSE, na.value = "grey80") + theme(legend.position = "none") + labs(title = NULL)
  
 #plot crop4
ImageDimPlot(meso.WS.Tumor, fov = "zoom4", group.by = "bin_Sarcomatoid", cols = c("gray80", yellow.hue), border.size = 0.1,
  border.color = "white", axes = FALSE, na.value = "grey80") + theme(legend.position = "none") + labs(title = NULL)
