library(Seurat)
library(pheatmap)
library(tidyverse)
library(cowplot)
library(dplyr)
set.seed(255)
setwd("~/Work/Collaborations/sp_2dpf_scRNAseq_analysis")
#1
sp_2dpf_1.data <- Read10X(data.dir = "./sp_2dpf_1")
sp_2dpf_1_barcodes <- read.delim(file = "./sp_2dpf_1/barcodes.tsv", stringsAsFactors = F, header = F)
sp_2dpf_1 <- CreateSeuratObject(sp_2dpf_1.data, project = "sp_2dpf_1", min.cells = 3, min.features = 300, meta.data = rownames(sp_2dpf_1_barcodes$v1)) #check values for cells and features

sp_2dpf_1a.data <- Read10X(data.dir = "./sp_2dpf_1a")
sp_2dpf_1a_barcodes <- read.delim(file = "./sp_2dpf_1a/barcodes.tsv", stringsAsFactors = F, header = F)
sp_2dpf_1a <- CreateSeuratObject(sp_2dpf_1.data, project = "sp_2dpf_1a", min.cells = 3, min.features = 350, meta.data = rownames(sp_2dpf_1a_barcodes$v1)) #check values for cells and features


#2
sp_2dpf_2.data <- Read10X(data.dir = "./sp_2dpf_2")
sp_2dpf_2_barcodes <- read.delim(file = "./sp_2dpf_2/barcodes.tsv", stringsAsFactors = F, header = F)
sp_2dpf_2 <- CreateSeuratObject(sp_2dpf_2.data, project = "sp_2dpf_2", min.cells = 3, min.features = 300, meta.data = rownames(sp_2dpf_2_barcodes$v1)) #check values for cells and features

sp_2dpf_2a.data <- Read10X(data.dir = "./sp_2dpf_2a")
sp_2dpf_2a_barcodes <- read.delim(file = "./sp_2dpf_2a/barcodes.tsv", stringsAsFactors = F, header = F)
sp_2dpf_2a <- CreateSeuratObject(sp_2dpf_2a.data, project = "sp_2dpf_2a", min.cells = 3, min.features = 500, meta.data = rownames(sp_2dpf_2a_barcodes$v1)) #check values for cells and features

#3
sp_2dpf_3.data <- Read10X(data.dir = "./sp_2dpf_3") #They call genes.tsv file barcode file for some reason, if its not present
sp_2dpf_3_barcodes <- read.delim(file = "./sp_2dpf_3/barcodes.tsv", stringsAsFactors = F, header = F)
sp_2dpf_3 <- CreateSeuratObject(sp_2dpf_3.data, project = "sp_2dpf_3", min.cells = 3, min.features = 200, meta.data = rownames(sp_2dpf_3_barcodes$v1)) #check values for cells and features



#Renamimg
sp_2dpf_1 <- RenameCells(sp_2dpf_1, add.cell.id = "sp_2dpf_1")
sp_2dpf_2 <- RenameCells(sp_2dpf_2, add.cell.id = "sp_2dpf_2")
sp_2dpf_3 <- RenameCells(sp_2dpf_3, add.cell.id = "sp_2dpf_3")
sp_2dpf_1a <-RenameCells(sp_2dpf_1a, add.cell.id = "sp_2dpf_1a")
sp_2dpf_2a <-RenameCells(sp_2dpf_2a, add.cell.id = "sp_2dpf_2a")


#Combination
sp_2dpf_list <- list(sp_2dpf_1, sp_2dpf_2, sp_2dpf_1a, sp_2dpf_3, sp_2dpf_2a)

for (i in 1:length(sp_2dpf_list)) {
  
  sp_2dpf_list[[i]] <- NormalizeData(sp_2dpf_list[[i]], 
                                  verbose = FALSE)
  sp_2dpf_list[[i]] <- FindVariableFeatures(sp_2dpf_list[[i]], 
                                         selection.method = "vst",
                                         nfeatures = 2000,
                                         verbose = FALSE)
  
}
sp_2dpf_anchors <- FindIntegrationAnchors(object.list = sp_2dpf_list, dims = 1:30)

sp_2dpf <- IntegrateData(anchorset = sp_2dpf_anchors, dims = 1:30)

#Clustering                                                                                         
sp_2dpf_integrated <- ScaleData(object = sp_2dpf)
sp_2dpf_integrated <- RunPCA(object = sp_2dpf_integrated, features = VariableFeatures(object = sp_2dpf_integrated))
sp_2dpf_integrated <- FindNeighbors(object = sp_2dpf_integrated, dims = 1:20)
sp_2dpf_integrated<- FindClusters(object = sp_2dpf_integrated, resolution = 1)

#UMAP
sp_2dpf_integrated_umap <- RunUMAP(object = sp_2dpf_integrated, dims = 1:20)
DimPlot(sp_2dpf_integrated_umap, group.by="orig.ident")
DimPlot(object = sp_2dpf_integrated_umap, reduction = "umap", label = TRUE)
DimPlot(sp_2dpf_integrated_umap, group.by="orig.ident", split.by = "orig.ident")

#Markers
markers <- FindAllMarkers(object = sp_2dpf_integrated_umap, assay = "RNA", only.pos = TRUE, min.pct = 0.01 )
genes <- read.delim("./names.tsv")
named_markers <- inner_join(markers, genes, by = "gene")
write.table(named_markers, file= "Sp_2dpf_marker_genes", quote= FALSE, sep = "\t")

#Extract all genes from dotplot
sp_2dpf_goi_info <- ggplot_build(DotPlot(sp_2dpf_integrated_umap, features = all))$plot$data
names <- read.delim("./names.tsv")
names <- names %>% rename("features.plot" = "gene", "name"="name")
sp_2dpf_goi_info <- inner_join(sp_2dpf_goi_info, names, by ="features.plot")
write.table(sp72_goi_info, file= "Sp_2dpf_all_genes", quote= FALSE, sep = "\t")

DotPlot(object = sp_2dpf_integrated_umap, features = "", col.min=0, scale.by = "size", cols = c("white", "red")) + coord_flip()

#Finding genes expressed in hindgut cells
sp_2dpf_hindgut_genes <- sp_2dpf_goi_info[sp_2dpf_goi_info$id == "Hindgut (1)" & sp_2dpf_goi_info$avg.exp.scaled >= 0.5,]

write.table(sp_2dpf_hindgut_genes, file= "sp_2dpf_hindgut_genes.tsv", quote= FALSE, sep = "\t")


#Extracting info from Sp-Pdx1 positive cells
pdx1 <- WhichCells(sp_2dpf_integrated_umap, idents = "Hindgut (1)", expression = `WHL22.169409` > 0)
sp_2dpf_integrated_umap <- AddModuleScore(sp_2dpf_integrated_umap, features = "WHL22.169409", seed = 255, name = "Pdx")
pdx1cells <- subset(sp_2dpf_integrated_umap, idents = "Hindgut (1)", subset = Pdx1 >= 1.5)
pdx1_all_genes_info <- as.data.frame(AverageExpression(pdx1cells, assays = "RNA", slot = "data", features = all_whl)) %>% rownames_to_column("whl")
pdx1_all_genes_info$cluster <-  "hindgut cells"
pdx1_expressed_genes_info <- pdx1_all_genes_info[pdx1_all_genes_info$all >= 0.5,]

write.table(pdx1_expressed_genes_info, file= "sp_2dpf_SpPdx_positive_hindgut_cells_genes.txt", quote= FALSE, sep = "\t")


