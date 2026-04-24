#load libraries
set.seed(1234)
library(cellchat2)
library(ggplot2)
library(tidyverse)

#Figure5B: cellchat comparative interaction weight bar plot for xenium meso TMA data.

#plot bar plot
compareInteractions(cellchat, show.legend = FALSE, group = c(1,2), measure = "weight")
