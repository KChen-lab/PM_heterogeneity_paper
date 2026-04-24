#load libraries
set.seed(1234)
library(cellchat2)
library(ggplot2)
library(tidyverse)

#Figure 5C: cellchat ligand-receptor pair plot for tumor cells (3) and CD8+ T cells (7).

#plot ligand-receptor dot plot
netVisual_bubble(cellchat, sources.use = c(3), targets.use = c(3,7), comparison = c(1,2), angle.x = 45)
