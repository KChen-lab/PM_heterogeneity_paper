#load libraries
set.seed(1234)
library(cellchat2)
library(ggplot2)
library(tidyverse)

#plot bar plot
compareInteractions(cellchat, show.legend = FALSE, group = c(1,2), measure = "weight")