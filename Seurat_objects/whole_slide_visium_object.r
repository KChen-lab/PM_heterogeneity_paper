#load libraries
set.seed(1234)
library(SeuratObject)
library(Seurat)
library(tidyverse)

#create seurat object from visium whole slide (WS)
vis.WS <- LoadData ("/path to data/",)