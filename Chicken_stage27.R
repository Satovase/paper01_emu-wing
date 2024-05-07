# R version 4.2.1
library(Seurat) # 4.3.0
library(patchwork) # ver 1.1.2
library(ggplot2) # 3.4.2
library(conflicted # 1.2.0
library(dplyr) # ver 1.1.2
library(cowplot) # ver 1.1.1

# Loading data
gg27sc.data <- Read10X(data.dir = "/Users/tanaka/run_count_GgHH27sc_GSM5319712/outs/filtered_feature_bc_matrix")
gg27sc <- CreateSeuratObject(counts=gg27sc.data, project="HH27sc", min.cells=3, min.features=200)

# QC
gg27sc[["percent.mt"]] <- PercentageFeatureSet(gg27sc, pattern = "MT-")
gg27sc <- subset(gg27sc, subset = percent.mt<15)
gg27sc <- subset(gg27sc, subset = nFeature_RNA < 1000 & nCount_RNA > 4000, invert = TRUE)
gg27sc <- subset(gg27sc, subset = nFeature_RNA < 300 & nCount_RNA > 1000, invert = TRUE)
gg27sc <- subset(gg27sc, subset = nFeature_RNA > 2000, invert = TRUE)
gg27sc <- subset(gg27sc, subset = nCount_RNA > 5000, invert = TRUE)

# Pre-processing

gg27sc <- NormalizeData(gg27sc, normalization.method="LogNormalize", scale.factor=10000)
gg27sc <- FindVariableFeatures(gg27sc, selection.method="vst", nfeatures=2000)
all.genes <- rownames(gg27sc)
gg27sc <- ScaleData(gg27sc, features=all.genes)
gg27sc <- RunPCA(gg27sc, features=VariableFeatures(object=gg27sc))

# TSNE
gg27sc <- FindNeighbors(gg27sc, dims=1:20)
gg27sc <- FindClusters(gg27sc, resolution=0.2)
gg27sc.t <- RunTSNE(gg27sc, dims=1:20)
DimPlot(gg27sc.t) # Fig. S13e

# Differential Expression Genes
all.markers <- FindAllMarkers(object = gg27sc.t)
write.csv(all.markers, file = "~/FindAllMarker_chicken27.csv") # Supplymentary table 5

# Dotplot
dotplot_genes <- c("TBX5", "PRRX1", "HAND2", "FN1", "TNNT3", "MYOD1", "MYOG", "MET", "PAX3", "LMO2", "SOX10")
p <- DotPlot(gg27sc.t, features = dotplot_genes,  cols = c("yellow", "red"), dot.scale = 8) +coord_flip() # Fig. S13f

# Muscle cluster analysis
muscle <- subset(gg27sc.t, ident = 3)
cell_barcords <- WhichCells(muscle)
gg27sc <- CreateSeuratObject(counts=gg27sc.data, project="HH27sc", min.cells=3, min.features=200)
subset <- subset(gg27sc, cells = cells_barcord)

# QC
subset <- NormalizeData(subset, normalization.method = "LogNormalize", scale.factor = 10000)
subset <- FindVariableFeatures(subset, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(subset)
subset <- ScaleData(subset, features = all.genes)
subset <- RunPCA(subset, eatures = VariableFeatures(object = subset))
ElbowPlot(subset, ndims = 20, reduction = "pca")
subset <- FindNeighbors(subset, dims = 1:10)
subset <- FindClusters(subset, resolution=1)
subset.u <- RunUMAP(subset, dims = 1:10)

# FeaturePlot
genes <-  c("PAX3", "HAND2", "PRRX1", "TBX5")
for (i in genes) {
  filename <- paste0(i, ".pdf")
  print(filename)
  p <- FeaturePlot(subset.u, features = i, cols = c("#e6e6e6", "blue"), pt.size = 1)
  pdf(file = filename)
  print(p)
  dev.off()
} # Fig. S13h
