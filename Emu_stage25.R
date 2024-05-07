# R version 4.2.1
library(Seurat) # 4.3.0
library(patchwork) # ver 1.1.2
library(ggplot2) # 3.4.2
library(conflicted # 1.2.0
library(dplyr) # ver 1.1.2
library(cowplot) # ver 1.1.1

# Loading data
emu25.data <- Read10X(data.dir = "~/run_count_EmuHH25sc_refseq/outs/filtered_feature_bc_matrix")
emu25 <- CreateSeuratObject(counts=emu25.data, project="EmuHH25_refseq", min.cells=3, min.features=200)
# 17050 features across 8710 samples

# QC
emu25[["percent.mt"]] <- PercentageFeatureSet(emu25, pattern="^KEF71-")

emu25 <- subset(emu25, subset = percent.mt<15)
emu25 <- subset(emu25, subset = nFeature_RNA > 200 & nCount_RNA > 3000)
emu25 <- subset(emu25, subset = nCount_RNA > 10000 & nFeature_RNA < 2000, invert = TRUE)
emu25 <- subset(emu25, subset = nCount_RNA > 45000 & nFeature_RNA < 7000, invert = TRUE)
emu25 <- subset(emu25, subset = nCount_RNA < 70000 & nFeature_RNA < 8500)
# 7363 samples

# Pre-processing
emu25 <- NormalizeData(emu25, normalization.method="LogNormalize", scale.factor=10000)
emu25 <- FindVariableFeatures(emu25, selection.method="vst", nfeatures=2000)
all.genes <- rownames(emu25)
emu25 <- ScaleData(emu25, features = all.genes)
emu25 <- RunPCA(emu25, features = VariableFeatures(object = emu25))

# TSNE
emu25.c <- FindNeighbors(emu25, dims=1:50)
emu25.c <- FindClusters(emu25.c, resolution=0.1)
emu25.t <- RunTSNE(emu25.c, dims = 1:50)
emu25.l <- RenameIdents(emu25.t, `0` = "LPM 1", `1` = "LPM 2", `2` = "Muscle",`3` = "Hematopietic cells", `4` = "Neural crest 1", `5` = "Neural crest 2")
DimPlot(emu25.l, label = TRUE) # Fig. 3d

# Differential Expression Genes
all.markers <- FindAllMarkers(object = emu25.t)
cluster <- all.markers[, "cluster"]
cluster <- gsub("2", "Muscle", cluster)
cluster <- gsub("1", "LPM 2", cluster)
cluster <- gsub("0", "LPM 1", cluster)
cluster <- gsub("3", "Hematopietic cells", cluster)
cluster <- gsub("4", "Neural crest 1", cluster)
cluster <- gsub("5", "Neural crest 2", cluster)
all.markers[, "cluster"] <- cluster
write.csv(all.markers, file = "~/FindAllMarker_emu25.csv") # Supplymentary table 2

# Dotplot
dotplot_genes <- c("TBX5", "PRRX1", "HAND2", "FN1", "TNNT3", "LBX1", "MYOD1", "MYOG", "MET", "PAX3", "LMO2", "CDH5", "SOX10")
DotPlot(emu25.l, features = dotplot_genes,  cols = c("yellow", "red"), dot.scale = 8) +coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust=1)) # Fig. 3e

# FeaturePlot
for (i in dotplot_genes) {
  filename <- paste0(i, ".pdf")
  print(filename)
  p <- FeaturePlot(emu25.t, features = i, cols = c("#e6e6e6", "blue"), pt.size = 1)
  pdf(file = filename)
  print(p)
  dev.off() # Fig. S11
}

# ViolinPlot
subset <- subset(emu25.t, ident = 2)
Pax3_expression = GetAssayData(object = subset, assay = "RNA", slot = "data")["PAX3",]
Hand2_expression = GetAssayData(object = subset, assay = "RNA", slot = "data")["HAND2",]

PH_pos_ids = names(which(Pax3_expression>0 & Hand2_expression>0))
H_pos_ids = names(which(Pax3_expression==0 & Hand2_expression>0))
P_pos_ids = names(which(Pax3_expression>0 & Hand2_expression==0))
PH_neg_ids = names(which(Pax3_expression==0 & Hand2_expression==0))

PH_pos_cells = subset(subset,cells=PH_pos_ids)
PH_pos_cells@meta.data[, "split_gene"] <- "Pax3+Hand2+cells" # 97 samples
P_cells = subset(subset,cells=P_pos_ids)
P_cells@meta.data[, "split_gene"] <- "Pax3+cells" # 328 samples
H_cells = subset(subset,cells=H_pos_ids)
H_cells@meta.data[, "split_gene"] <- "Hand2+cells" # 11 samples
PH_neg_cells = subset(subset,cells=PH_neg_ids)
PH_neg_cells@meta.data[, "split_gene"] <- "Pax3-_Hand2-" # 69 samples

genes <- c("HAND2", "PAX3", "TNNT3", "MYOD1", "TBX5", "MET", "LBX1", "CASP10", "APAF1", "BAK1")
split_cells <- merge(P_cells, PH_pos_cells)
order <- c("Pax3+cells", "Pax3+Hand2+cells")
split_cells@meta.data$split_gene <- factor(x = split_cells@meta.data$split_gene, levels = order)
for (i in genes) {
  filename <- paste0(i, "_musclecluster_Pax3-Hand2.pdf")
  print(filename)
  p <- VlnPlot(split_cells, features = i, group.by = "split_gene", cols = c("#FFFF00", "#FFA500", "#FF00FF"))
  pdf(file = filename)
  print(p)
  dev.off() # Fig. 3g, 5a
}

# Statistical processing
Idents(object = split_cells) <- split_cells@meta.data$'split_gene'
markers <- FindMarkers(split_cells, ident.1 = "Pax3+cells", ident.2 = "Pax3+Hand2+cells", logfc.threshold = 0)
head(markers)
write.csv(markers, file = "~/FindMarker_emu25.csv") # Supplymentary table 3

# Muscle cluster analysis
muscle <- subset(emu25.t, ident = 2)
cell_barcords <- WhichCells(muscle)
emu25 <- CreateSeuratObject(counts=gg24sc.data, project="HH24sc", min.cells=3, min.features=200)
subset <- subset(emu25, cells = cells_barcord)

# QC
subset <- NormalizeData(subset, normalization.method = "LogNormalize", scale.factor = 10000)
subset <- FindVariableFeatures(subset, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(subset)
subset <- ScaleData(subset, features = all.genes)
subset <- RunPCA(subset, eatures = VariableFeatures(object = subset))
ElbowPlot(subset, ndims = 20, reduction = "pca")
subset <- FindNeighbors(subset, dims = 1:10)
subset <- FindClusters(subset, resolution=0.2)
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
} # Fig. S12


