#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(viridis)

#color gradient
plasma_grad <- plasma(40)

#plot with correct segmentation boundaries
p1 <- ImageDimPlot(meso, fov = "fov", group.by = "refined_tumor_cluster", cols = cell.colors, border.size = 0.1) +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())

#plot with incorrect segmentation boundaries
p2 <- ImageFeaturePlot(meso, fov = "fov", features = "pseudotime", cols = plasma_grad, dark.background = TRUE, 
                 border.size = 0, boundaries = "segmentations") +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())

#add p1 (x,y) coordinates to p2
p2[[1]]$data <- p2[[1]]$data %>%
  left_join(p1[[1]]$data %>% select(cell, x_correct = x, y_correct = y),
            by = "cell") %>%
  mutate(x = x_correct, y = y_correct,
         x_correct = NULL, y_correct = NULL)
