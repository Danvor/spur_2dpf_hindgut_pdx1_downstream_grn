Sys.setenv(LANG = "en")
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(plotrix)
library(RColorBrewer)
set.seed(255)
setwd("~/Work/Collaborations/sp_2dpf_umap")
DimPlot(object = sp_2dpf_umap, reduction = "umap", label = FALSE, repel = TRUE) #+ theme(legend.position="none")

#WHL_names.tsv file linking every WHL id to gene name

#Finding genes expressed in hindgut cells
sp_2dpf_hindgut_genes <- sp_2dpf_goi_info[sp_2dpf_goi_info$id == "Hindgut (1)" & sp_2dpf_goi_info$avg.exp.scaled >= 0.5,]

write.table(sp_2dpf_hindgut_genes, file= "sp_2dpf_hindgut_genes.tsv", quote= FALSE, sep = "\t")


#Extracting info from Sp-Pdx1 positive cells
pdx1 <- WhichCells(sp_2dpf_umap, idents = "Hindgut (1)", expression = `WHL22.169409` > 0)


sp_2dpf_umap <- AddModuleScore(sp_2dpf_umap, features = "WHL22.169409", seed = 255, name = "Pdx")
FeaturePlot(pdx1cells, features = "Pdx1", label = TRUE, repel = TRUE, pt.size = 1, order = TRUE) +
  scale_colour_gradientn(colors = c("White","Red","DarkGreen","Purple"))

pdx1cells <- subset(sp_2dpf_umap, idents = "Hindgut (1)", subset = Pdx1 >= 1.5)
pdx1_all_genes_info <- as.data.frame(AverageExpression(pdx1cells, assays = "RNA", slot = "data", features = all_whl)) %>% rownames_to_column("whl")
pdx1_all_genes_info$cluster <-  "hindgut cells"
pdx1_expressed_genes_info <- pdx1_all_genes_info[pdx1_all_genes_info$all >= 0.5,]
write.table(pdx1_expressed_genes_info, file= "sp_2dpf_SpPdx_positive_hindgut_cells_genes.txt", quote= FALSE, sep = "\t")
