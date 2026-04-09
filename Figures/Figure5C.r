#load libraries
set.seed(1234)
library(cellchat2)
library(ggplot2)
library(tidyverse)

#plot ligand-receptor dot plot
netVisual_bubble(cellchat, sources.use = c(6), targets.use = c(3,6), comparison = c(1,2), angle.x = 45)
