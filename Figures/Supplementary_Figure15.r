#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#panel A
SpatialFeaturePlot(vis.WS, features = c("KEGG_RETINOL_METABOLISM"))

#panel B
SpatialFeaturePlot(vis.WS, features = c("KEGG_PENTOSE_PHOSPHATE_PATHWAY"))

#panel C
SpatialFeaturePlot(vis.WS, features = c("KEGG_PURINE_METABOLISM"))

#panelD
SpatialFeaturePlot(vis.WS, features = c("KEGG_NICOTINATE_AND_NICOTINAMIDE_METABOLISM"))
