#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#cell cycle genes 
s.genes.xen <- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM7", "MCM4", "MCM6", "RRM1", "UNG", "UHRF1",
                 "HELLS", "MSH2", "RAD51", "RRM2", "EXO1", "BLM", "CASP8AP2")
g2m.genes.xen <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "CKS2", "MKI67",
                   "TMPO", "CENPF", "PIMREG", "SMC4", "CCNB2", "CKAP2", "AURKB", "AURKA", "BUB1",
                   "KIF11", "ANP32E", "CDCA3", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
                   "CDCA2", "HMMR", "PSRC1", "CKAP5", "CTCF", "NEK2", "CBX5", "CENPA")

#transition signature
trans.sig <- c("SNAI2", "TWIST1", "TWIST2", "MYBL2")

#assign cell cycle scores
meso.WS <- CellCycleScoring(meso.WS, s.features = s.genes.xen, g2m.features = g2m.genes.xen, set.ident = TRUE)

#score transition signaturr with module score
meso.WS <- AddModuleScore(meso.WS, features=list(trans.sig), name= "Transition_signature", ctrl = 50)

#statistics
g2m.cor <- cor.test(meso.WS$Transition_signature1, meso.WS$G2M.Score, method = "spearman", exact = FALSE)
g2m.cor$estimate
g2m.cor$p.value

g2m.rho <- round(unname(g2m.cor$estimate), 3)
g2m.p.label <- ifelse(g2m.cor$p.value < 2.2e-16,
  "p < 2.2e-16", paste0("p = ", signif(g2m.cor$p.value, 3)))
g2m.label <- paste0("rho = ", rho, "\n", g2m.p.label)

s.cor <- cor.test(meso.WS$Transition_signature1, meso$S.Score, method = "spearman")
s.cor$estimate
s.cor$p.value

s.rho <- round(unname(s.cor$estimate), 3)
s.p.label <- ifelse(s.cor$p.value < 2.2e-16,
                      "p < 2.2e-16", paste0("p = ", signif(s.cor$p.value, 3)))
s.label <- paste0("rho = ", rho, "\n", s.p.label)


#cell cycle correlated with transition
g1 <- FeatureScatter(meso.WS, feature1 = "Transition_signature1", feature2 = "G2M.Score", group.by = "Phase") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = g2m.label)

g2 <- FeatureScatter(meso.WS, feature1 = "Transition_signature1", feature2 = "S.Score", group.by = "Phase") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = s.label)
g1 + g2
