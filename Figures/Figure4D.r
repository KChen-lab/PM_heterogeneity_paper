#load libraries
set.seed(1234)
library(Seurat)
library(ggplot2)
library(tidyverse)

#create tumor only subcluster
Idents(meso) <- "refined_tumor_cluster"
mesoTumor <- subset(meso, idents = c("Epithelioid tumor", "Epithelioid proliferating tumor",
                                     "Sarcomatoid proliferating tumor", "Sarcomatoid tumor"))

#create df of tumor state counts
df <- mesoTumor@meta.data %>%
  filter(!is.na(patient_number)) %>%
  group_by(patient_number, refined_tumor_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patient_number) %>%
  mutate(prop = n / sum(n))

#plot is
tumor.colors <- c(
  "Sarcomatoid tumor" = "#006896", "Epithelioid tumor" = "#e84359", "Sarcomatoid proliferating tumor" = "#5abacc", 
  "Epithelioid proliferating tumor" = "#ff8fb3")

#plot bar plot
ggplot(df, aes(x = patient_number, y = prop, fill = refined_tumor_cluster)) + 
  geom_bar(stat = "identity", width = 0.8) + coord_flip() +
  scale_fill_manual(values = tumor.colors) +
  labs(x = "Sample", y = "Proportion", fill = "Tumor State") + theme_cowplot()
