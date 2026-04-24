#load libraries
library(EnhancedVolcano)

#Figure 5F: Volcano plot for significantly expressed genes between epithelioid and biphasic TME ROIs from GeoMx DSP WTA data.

#genes to label
genes_to_label <- c("CD27", "ITGB4", "HAVCR2", "CD8A",
                    "CXCL11", "LAG3", "PRF1",
                    "FOS", "JUN", "SPARC", "CCN2", "GZMB", "CXCL11",
                    "NKG7")

#make genes italic
lab.italics <- paste0("italic('", rownames(histoTMECombined_topGenes), "')")
selectLab.italics <- paste0("italic('", genes_to_label, "')")

#create key values for volcano plot
keyvals <- ifelse((histoTMECombined_topGenes$logFC < -1) & (histoTMECombined_topGenes$adj.P.Val < 0.05), '#6495ED',
                  ifelse((histoTMECombined_topGenes$logFC > 1) & (histoTMECombined_topGenes$adj.P.Val < 0.05), '#B22222',
                         ifelse((histoTMECombined_topGenes$logFC > 1) & (histoTMECombined_topGenes$adj.P.Val > 0.05), 'lightgray',
                                ifelse((histoTMECombined_topGenes$logFC < -1) & (histoTMECombined_topGenes$adj.P.Val > 0.05), 'lightgray',
                                       ifelse((abs(histoTMECombined_topGenes$logFC) < 1) & (histoTMECombined_topGenes$adj.P.Val < 0.05), 'gray',
                                              'darkgray')))))

names(keyvals)[keyvals == "#B22222"] <- "p.adj < 0.05 & logFC > 1"
names(keyvals)[keyvals == "#6495ED"] <- "p.adj < 0.05 & logFC < -1"
names(keyvals)[keyvals == "darkgray"] <- "Not Significant"
names(keyvals)[keyvals == "lightgray"] <- "Log2FC only"
names(keyvals)[keyvals == "gray"] <- "Adjusted p-value"

#plot volcano
EnhancedVolcano(histoTMECombined_topGenes, lab = lab.italics, x = "logFC", y = "adj.P.Val",
                selectLab = selectLab.italics, colAlpha = 1, parseLabels = TRUE,
                colCustom = keyvals, legendPosition = "right",
                FCcutoff = 1, pCutoff = 0.05, labSize = 5,
                xlim = c(-3.1, 3), ylim = c(0, 13),
                subtitle = "Adjusted p-value < 0.05, Absolute LogFC >1",
                ylab = bquote(~-Log[10]~Adjusted~italic(P)),
                xlab = bquote(~Log~ "FC"),
                title = "ePM Tumor vs. bPM Tumor",
                drawConnectors = TRUE, widthConnectors = 0.5, max.overlaps = Inf)
