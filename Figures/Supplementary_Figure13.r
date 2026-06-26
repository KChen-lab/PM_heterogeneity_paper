#load libraries
set.seed(1234)
library(CellChat2)
library(ggplot2)
library(tidyverse)

#plot stacked bar plot
rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
