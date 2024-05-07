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
# pca 結果書き出し
pdf(file = "/Users/tanaka/Analysis_GgHH24sc/published_data/pca/umap.pdf")
print(DimPlot(gg24sc.u))
dev.off()
pdf(file = "/Users/tanaka/Analysis_GgHH24sc/published_data/pca/tsne.pdf")
print(DimPlot(gg24sc.t))
dev.off()

# UMI 結果書き出し
pdf(file = "/Users/tanaka/Analysis_GgHH24sc/published_data/pca/umap-umi.pdf")
print(FeaturePlot(gg24sc.u, features = "nCount_RNA", cols = c("#009BBF", "red"), pt.size = 1))
dev.off()
pdf(file = "/Users/tanaka/Analysis_GgHH24sc/published_data/pca/tsne-umi.pdf")
print(FeaturePlot(gg24sc.t, features = "nCount_RNA", cols = c("#009BBF", "red"), pt.size = 1))
dev.off()

# featureplot書き出しルーティーン
# MYOGなし
genes <-  c("HAND2", "PAX3", "PAX7", "PRRX1", "MET", "MYOD1", "TNNT3", "EPHA3", "CD44", "TBX5")
pcluster_genes <- c("SOX10", "CDH5", "LMO2", "NKX2-5", "MET", "MYOG", "TNNT3", "COL4A2", "FN1", "ALX4", "PTCH1", "HAND2", "PRRX1")
pfeature_genes <- c("PRRX1", "HAND2", "PTCH1", "HOXD10", "MEIS1", "ALX4", "FN1", "COL4A2", "MYOG", "LBX1", "MYOD1", "LMO2", "CDH5", "SOX10")
LPM_genes <- c("PRRX1", "HAND1", "HAND2", "GATA4", "HOXB6", "BMP4", "TBX5", "FOXF1")

setwd("/Users/tanaka/Analysis_GgHH24sc/published_data/feat/u")
for (i in genes) {
  filename <- paste0(i, ".pdf")
  print(filename)
  p <- FeaturePlot(gg24sc.u, features = i, cols = c("#e6e6e6", "blue"), pt.size = 1)
  pdf(file = filename)
  print(p)
  dev.off()
}
setwd("/Users/tanaka/Analysis_GgHH24sc/published_data/feat/t")
for (i in genes) {
  filename <- paste0(i, ".pdf")
  print(filename)
  p <- FeaturePlot(gg24sc.t, features = i, cols = c("#e6e6e6", "blue"), pt.size = 1)
  pdf(file = filename)
  print(p)
  dev.off()
}

## Dotplot
dotplot_genes <- c("TBX5", "PRRX1", "HAND2", "FN1", "TNNT3", "LBX1", "MYOD1", "MYOG", "MET", "PAX3", "LMO2", "CDH5", "SOX10")
p <- DotPlot(gg24sc.t, features = dotplot_genes,  cols = c("yellow", "red"), dot.scale = 8) +coord_flip()
pdf(file = "/Users/tanaka/Analysis_GgHH24sc/published_data/dotplot/dotplot.pdf")
print(p)
dev.off()

## ViolinPlot
HH24_muscle <- subset(gg24sc.t, ident = 5)
DimPlot(HH24_muscle)
cell_barcords <- WhichCells(HH24_muscle)
saveRDS(cell_barcords, file = "/Users/tanaka/Analysis_GgHH24sc/saveRDS/240326_muscle_subset_cellbarcord")

subset <- HH24_muscle
genes <- c("HAND2", "PAX3", "TNNT3", "MYOD1", "MYOG", "MET", "LBX1")
setwd("/Users/tanaka/Analysis_GgHH24sc/published_data/feat/muscle_cluster")
for (i in genes) {
  filename <- paste0(i, ".pdf")
  print(filename)
  p <- FeaturePlot(subset, features = i, cols = c("#e6e6e6", "blue"), pt.size = 1)
  pdf(file = filename)
  print(p)
  dev.off()
}

Pax3_expression = GetAssayData(object = subset, assay = "RNA", slot = "data")["PAX3",]
Hand2_expression = GetAssayData(object = subset, assay = "RNA", slot = "data")["HAND2",]

PH_pos_ids = names(which(Pax3_expression>0 & Hand2_expression>0))
H_pos_ids = names(which(Pax3_expression==0 & Hand2_expression>0))
P_pos_ids = names(which(Pax3_expression>0 & Hand2_expression==0))
PH_neg_ids = names(which(Pax3_expression==0 & Hand2_expression==0))
Pax3_ids = names(which(Pax3_expression>0))
Hand2_ids = names(which(Hand2_expression>0))

PH_pos_cells = subset(subset,cells=PH_pos_ids)
# No cells found
P_pos_cells = subset(subset,cells=P_pos_ids)
# 50 samples
H_pos_cells = subset(subset,cells=H_pos_ids)
# No cells found
PH_neg_cells = subset(subset,cells=PH_neg_ids)
# 96 samples
Pax3_pos_cells = subset(subset,cells=Pax3_ids)
# 50 samples
Hand2_pos_cells = subset(subset,cells=Hand2_ids)
# オブジェクト 'Hand2_pos_cells' がありません 
