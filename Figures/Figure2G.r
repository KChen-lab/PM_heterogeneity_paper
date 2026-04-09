#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#colors
module.colors <- c(
  "Epithelioid_signature1" = "#e84359", "Transition_signature1" = "#a29bfe",
  "Sarcomatoid_signature1" = "#006896")
  
  p1 <- ggplot(
  morpho_long %>% filter(metric == "circularity"),
  aes(x = bin, y = value, fill = bin)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.75, aes(color = bin), linewidth = 0.9) +
  scale_fill_manual(values = c("#e84359", "#a29bfe", "#006896")) +
  scale_color_manual(values = c("#d7304f", "#8e8df0", "#005f73")) +
  labs(x = "Pseudotime bin", y = "Circularity") +
  theme(legend.position = "right") + theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  stat_summary(fun = median, geom = "point", size = 2, aes(color = bin)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1)))

#statisics
kruskal.test(circularity ~ bin, data = meta)
pairwise.wilcox.test(meta$circularity, meta$bin, p.adjust.method = "bonferroni")

#plot significance
p1 + stat_compare_means(method = "kruskal.test", label.y = 1.07, size = 3) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("Epithelioid", "Transition"),
      c("Transition", "Sarcomatoid"),
      c("Epithelioid", "Sarcomatoid")),
    label = "p.signif",
    step.increase = 0.08)
