#load libraries
set.seed(1234)
library(fgsea)
library(ggplot2)
library(tidyverse)

#dotplot function
fgseaDotPlot <- function(fgsea_results,
                         condition1 = "Upregulated",
                         condition2 = "Downregulated",
                         title = "Pathway Enrichment by Condition") {
  
  #add a sign column to splits pathways by enrichment to annotate by condition
  fgsea_results <- fgsea_results %>%
    mutate(.sign = ifelse(NES > 0, condition1, condition2),
           log10padj = -log10(padj))
  
  #filter out rows with NA 
  fgsea_results <- fgsea_results %>%
    filter(!is.na(NES), !is.na(pval), !is.na(padj))
  
  #get top pathways for each condition
  topPathways <- fgsea_results %>%
    group_by(.sign) %>%
    slice_min(order_by = padj, n = 11) %>%
    ungroup() %>%
    mutate(pathway = factor(pathway, levels = rev(unique(pathway))))
  
  #plot the top enriched pathways for each condition
  ggplot(topPathways, aes(x = NES, y = pathway, size = size, color = padj)) +
    geom_point() +
    facet_grid(. ~ .sign, scales = "free_y", space = "free_y") +
    scale_color_gradient(low = "#6495ED", high = "#B22222") +
    theme_light() +
    labs(title = title,
         x = "Normalized Enrichment Score (NES)",
         y = "Pathways",
         size = "Number of Genes\n in Pathway List",
         color = "p-adjusted value")
}

#plot dot plot
fgseaDotPlot(fgseaHistoTumorCom, "Epithelioid", "Biphasic", "Pathway Enrichment in Epithelioid vs. Biphasic Tumors")