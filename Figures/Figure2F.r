#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#figure 2F: area plot for module scores of tumor state gene signatures.
#source data (smooth.df) is in the supplementary data of the paper.

#colors
module.colors <- c(
  "Epithelioid_signature1" = "#e84359", "Transition_signature1" = "#a29bfe",
  "Sarcomatoid_signature1" = "#006896")
  
#plot area plot
ggplot(smooth.df, aes(x = pseudotime, y = score, fill = state)) +
  geom_area(alpha = 0.75, position = "identity") +
  theme_cowplot() +
  scale_fill_manual(values = module.colors) +
  labs(x = "Pseudotime", y = "Module Score", fill = "Tumor state") +
  theme(legend.position = "none")
