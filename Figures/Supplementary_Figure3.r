#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#plot VIM expression
SpatialFeaturePlot(meso, features = c("VIM"))
