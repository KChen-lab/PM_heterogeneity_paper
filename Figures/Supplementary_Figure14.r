#load libraries
set.seed(1234)
library(pheatmap)
library(ggplot2)
library(tidyverse)

#read z scores in
heat.mat.z <- read_csv("visium_whole_slide_tumor_metab_zscore.csv")
heat.mat.z <- as.matrix(heat.mat.z)

#plot heatmap
pheatmap(
  heat.mat.z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 7)

