# R version 4.2.1
library(Seurat) # 4.3.0
library(patchwork) # ver 1.1.2
library(ggplot2) # 3.4.2
library(conflicted # 1.2.0
library(dplyr) # ver 1.1.2
library(cowplot) # ver 1.1.1

# Loading data
gg24sc.data <- Read10X(data.dir = "~/run_count_GgHH24sc_GSM5319711/outs/filtered_feature_bc_matrix")
gg24sc <- CreateSeuratObject(counts=gg24sc.data, project="HH24sc", min.cells=3, min.features=200)

# QC
gg24sc[["percent.mt"]] <- PercentageFeatureSet(gg24sc, pattern = "MT-")
gg24sc <- subset(gg24sc, subset = percent.mt<15)
gg24sc <- subset(gg24sc, subset = nFeature_RNA < 2000 & nCount_RNA > 9000, invert = TRUE)
gg24sc <- subset(gg24sc, subset = nFeature_RNA < 1300 & nCount_RNA > 5000, invert = TRUE)
gg24sc <- subset(gg24sc, subset = nFeature_RNA < 950 & nCount_RNA > 3000, invert = TRUE)
gg24sc <- subset(gg24sc, subset = nFeature_RNA < 700 & nCount_RNA > 2100, invert = TRUE)
gg24sc <- subset(gg24sc, subset = nFeature_RNA > 2500, invert = TRUE)
gg24sc <- subset(gg24sc, subset = nCount_RNA > 10000, invert = TRUE)

# Pre-processing
gg24sc <- NormalizeData(gg24sc, normalization.method="LogNormalize", scale.factor=10000)
gg24sc <- FindVariableFeatures(gg24sc, selection.method="vst", nfeatures=2000)
all.genes <- rownames(gg24sc)
gg24sc <- ScaleData(gg24sc, features=all.genes)
gg24sc <- RunPCA(gg24sc, features=VariableFeatures(object=gg24sc))

# TSNE
gg24sc <- FindNeighbors(gg24sc, dims=1:20)
gg24sc <- FindClusters(gg24sc, resolution=0.2)
gg24sc.t <- RunTSNE(gg24sc, dims=1:20)
DimPlot(gg24sc.t) # Fig. S13a

# Differential Expression Genes
all.markers <- FindAllMarkers(object = gg24sc.t)
write.csv(all.markers, file = "~/FindAllMarker_chicken24.csv") # Supplymentary table 4

# Dotplot
dotplot_genes <- c("TBX5", "PRRX1", "HAND2", "FN1", "TNNT3", "LBX1", "MYOD1", "MYOG", "MET", "PAX3", "LMO2", "CDH5", "SOX10")
DotPlot(gg24sc.t, features = dotplot_genes,  cols = c("yellow", "red"), dot.scale = 8) +coord_flip() # Fig. S13b

# Muscle cluster analysis
muscle <- subset(gg24sc.t, ident = 5)
cell_barcords <- WhichCells(muscle)
gg24sc <- CreateSeuratObject(counts=gg24sc.data, project="HH24sc", min.cells=3, min.features=200)
subset <- subset(gg24sc, cells = cells_barcord)

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
} # Fig. S13d
