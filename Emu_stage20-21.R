# R version 4.2.1
library(Seurat) # 4.3.0
library(patchwork) # ver 1.1.2
library(ggplot2) # 3.4.2
library(conflicted # 1.2.0
library(dplyr) # ver 1.1.2
library(cowplot) # ver 1.1.1

# Loading data
emu2021.data <- Read10X(data.dir = "~/filtered_feature_bc_matrix")
emu2021 <- CreateSeuratObject(counts=emu2021.data, project="EmuHH2021_refseq", min.cells=3, min.features=200)
# 16868 features across 3871 samples

# QC
emu2021[["percent.mt"]] <- PercentageFeatureSet(emu2021, pattern="^KEF71-")
emu2021 <- subset(emu2021, subset = percent.mt<15)
emu2021 <- subset(emu2021, subset = nFeature_RNA > 200 & nCount_RNA > 5000)
emu2021 <- subset(emu2021, subset = nFeature_RNA < 1000, invert = TRUE)
emu2021 <- subset(emu2021, subset = nCount_RNA > 12000 & nFeature_RNA < 2500, invert = TRUE)
emu2021 <- subset(emu2021, subset = nCount_RNA > 20000 & nFeature_RNA < 3000, invert = TRUE)
emu2021 <- subset(emu2021, subset = nCount_RNA > 45000 & nFeature_RNA < 7000, invert = TRUE)
emu2021 <- subset(emu2021, subset = nCount_RNA < 70000 & nFeature_RNA < 8500)
# 2696 samples

# Pre-processing
emu2021 <- NormalizeData(emu2021, normalization.method="LogNormalize", scale.factor=10000)
emu2021 <- FindVariableFeatures(emu2021, selection.method="vst", nfeatures=2000)
all.genes <- rownames(emu2021)
emu2021 <- ScaleData(emu2021, features=all.genes)
emu2021 <- RunPCA(emu2021, features=VariableFeatures(object=emu2021))

# TSNE
emu2021 <- FindNeighbors(emu2021, dims=1:50)
emu2021 <- FindClusters(emu2021, resolution=0.2)
emu2021.t <- RunTSNE(emu2021, dims = 1:50)
emu2021.l <- RenameIdents(emu2021.t, `0` = "LPM 1", `1` = "LPM 2", `2` = "Muscle",`3` = "LPM 3", `4` = "Neural crest", `5` = "Intermerdiate mesoderm", `6` = "Neuron", `7` = "Hematopietic cells", `8` = "Ectoderm", `9` = "Notochord")
order <- c("Ectoderm", "Notochord", "Neuron", "Neural crest", "Hematopietic cells", "Intermerdiate mesoderm", "Muscle", "LPM 3", "LPM 2", "LPM 1")
DimPlot(emu2021.l, order = order) # Fig.3a
        
# Differential Expression Genes
all.markers <- FindAllMarkers(object = emu2021.t)
cluster <- all.markers[, "cluster"]
cluster <- gsub("2", "Muscle", cluster)
cluster <- gsub("1", "LPM 2", cluster)
cluster <- gsub("0", "LPM 1", cluster)
cluster <- gsub("3", "LPM 3", cluster)
cluster <- gsub("4", "Neural crest", cluster)
cluster <- gsub("5", "Intermerdiate mesoderm", cluster)
cluster <- gsub("6", "Neuron", cluster)
cluster <- gsub("7", "Hematopietic cells", cluster)
cluster <- gsub("8", "Ectoderm", cluster)
cluster <- gsub("9", "Notochord", cluster)
all.markers[, "cluster"] <- cluster
write.csv(all.markers, file = "~/FindAllMarker_emu2021.csv") # supplymentaxry table 1

# Dotplot
dotplot_genes <- c("TWIST1", "EBF3", "TBX5", "PRRX1", "HAND2", "LBX1", "MYOD1", "PAX3", "TNNT3", "MYOG", "MET", "OSR1", "PAX2", "LMO2", "CDH5", "SOX10", "ELAVL4", "SHH", "TBXT", "WNT6")        
levels(emu2021.l) <- c("LPM 1", "LPM 2", "LPM 3", "Muscle", "Intermerdiate mesoderm", "Hematopietic cells", "Neural crest", "Neuron", "Notochord", "Ectoderm")
DotPlot(emu2021.l, features = dotplot_genes,  cols = c("yellow", "red"), dot.scale = 8) +coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust=1)) # Fig. 3b

# FeaturePlot
genes <- c("TWIST1", "EBF3", "TBX5", "PRRX1", "HAND2", "OSR1", "PAX2", "LBX1", "MYOD1", "PAX3", "TNNT3", "MYOG", "MET", "SOX10", "ELAVL4", "LMO2", "CDH5", "WNT6", "SHH", "TBXT")
for (i in genes) {
  filename <- paste0(i, ".pdf")
  print(filename)
  p <- FeaturePlot(emu2021.t, features = i, cols = c("#e6e6e6", "blue"), pt.size = 1)
  pdf(file = filename)
  print(p)
  dev.off() # Fig. S9
}

subset <- subset(emu2021.t, ident = 2)
genes <- c("HAND2", "PAX3", "TNNT3", "MYOD1", "MYOG", "MET", "LBX1")
for (i in genes) {
  filename <- paste0(i, ".pdf")
  print(filename)
  p <- FeaturePlot(subset, features = i, cols = c("#e6e6e6", "blue"), pt.size = 1)
  pdf(file = filename)
  print(p)
  dev.off() # Fig. S10
}



        
