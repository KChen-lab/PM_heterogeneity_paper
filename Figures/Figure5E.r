#load libraries
library(EnhancedVolcano)

#genes to label
genes_to_label <- c("MSLN", "CALB2", "WT1", "CDH3", "AXL", "MYBL2",
                    "CD70", "GAS6", "LTBP2", "GAS6", "ITGB4",
                    "CD8A", "LAG3", "S100A4", "PRF1", "KDR", "KRT19",
                    "FOS", "JUN", "SPARC", "CCN2", "GZMB",
                    "CLDN15", "COBLL1", "RERG", "CCL8",
                    "PARP14", "IFIT3")

#make genes italic
lab.italics <- paste0("italic('", rownames(histoTumorCombined_topGenes), "')")
selectLab.italics <- paste0("italic('", genes_to_label, "')")

#create key values for volcano plot
keyvals <- ifelse((histoTumorCombined_topGenes$logFC < -1) & (histoTumorCombined_topGenes$adj.P.Val < 0.05), '#6495ED',
                  ifelse((histoTumorCombined_topGenes$logFC > 1) & (histoTumorCombined_topGenes$adj.P.Val < 0.05), '#B22222',
                         ifelse((histoTumorCombined_topGenes$logFC > 1) & (histoTumorCombined_topGenes$adj.P.Val > 0.05), 'lightgray',
                                ifelse((histoTumorCombined_topGenes$logFC < -1) & (histoTumorCombined_topGenes$adj.P.Val > 0.05), 'lightgray',
                                       ifelse((abs(histoTumorCombined_topGenes$logFC) < 1) & (histoTumorCombined_topGenes$adj.P.Val < 0.05), 'gray',
                                              'darkgray')))))

names(keyvals)[keyvals == "#B22222"] <- "p.adj < 0.05 & logFC > 1"
names(keyvals)[keyvals == "#6495ED"] <- "p.adj < 0.05 & logFC < -1"
names(keyvals)[keyvals == "darkgray"] <- "Not Significant"
names(keyvals)[keyvals == "lightgray"] <- "Log2FC only"
names(keyvals)[keyvals == "gray"] <- "Adjusted p-value"

#plot volcano
EnhancedVolcano(histoTumorCombined_topGenes, lab = lab.italics, x = "logFC", y = "adj.P.Val",
                selectLab = selectLab.italics, colAlpha = 1, parseLabels = TRUE,
                colCustom = keyvals, legendPosition = "right",
                FCcutoff = 1, pCutoff = 0.05, labSize = 5,
                xlim = c(-3.1, 3), ylim = c(0, 13),
                subtitle = "Adjusted p-value < 0.05, Absolute LogFC >1",
                ylab = bquote(~-Log[10]~Adjusted~italic(P)),
                xlab = bquote(~Log~ "FC"),
                title = "ePM Tumor vs. bPM Tumor",
                drawConnectors = TRUE, widthConnectors = 0.5, max.overlaps = Inf)

