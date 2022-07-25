#!/usr/bin/env Rscript
# module load R/4.0.2

# Load libraries:
library(melange)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(clustree)
library(stringr)

# Import RSEM counts using Melange:
cell.names <- scan("ListOfCellsWithRSEMgeneresults", character(), quote = "")
cell.rsem.file.locations <- scan("ListOfRSEMgeneresults", character(), quote = "")

# Create a Melange object containing collated quantification data:
rsem.melange <- LoadRSEM(cell.names, cell.rsem.file.locations)

# Save Melange object for future use:
saveRDS(rsem.melange, file = "rsem_melange.rds")

# Load saved Melange object to R workspace:
readRDS(file = "rsem_melange.rds")

# NOTE
# Authors of the reference publication state that cells with less than 900 genes detected were excluded for quality control filtering.
# In addition, genes that were not detected in at least 3 of the cells after this filtering were also removed from further analysis.
# For protein-coding genes we will keep this limit. For lncRNAs


#  --- Edit --- Create a new Seurat object using the data from melange:
# with protein-coding genes:
seurat.object <- CreateSeuratObject(counts = rsem.melange@data, project = "MySeuratProject", assay = "RNA", min.cells = 3, min.features = 900)
ncol(seurat.object)

# without protein-coding genes:
#seurat.object <- CreateSeuratObject(counts = rsem.melange@data, project = "MySeuratProject", assay = "RNA", min.cells = 3, min.features = 300)
#ncol(seurat.object)

# Write the number of genes in a table:
write.table(seurat.object$nFeature_RNA, file="SeuratGenesPerCell", quote=FALSE, sep='\t', col.names = FALSE)

# Convert a Melange object directly into a Seurat assay (currently not used):
#seurat.assay <- SeuratAssayFromMelange(rsem.melange)

# Identify mitochondrial reads:
seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, pattern = "^MT-")
seurat.object[["percent.mt"]][is.na(seurat.object[["percent.mt"]])] <- 0

# Alternatively:
# nFeature_RNA is the number of genes detected in each cell (d=200-2500).
# nCount_RNA is the total number of molecules detected within a cell.
#seurat.object <- subset(seurat.object, subset = nFeature_RNA > 600 & nFeature_RNA < 7000 & percent.mt < 5)
#ncol(seurat.object)

# Plot main features as violin plots:
tiff("ViolinPlots_Features.tiff", units="in", width=5, height=5, res=300)
VlnPlot(seurat.object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

# Plot main features as scatter plots:
plot1 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

tiff("ScatterPlot_FilteredFeatures.tiff", units="in", width=5, height=5, res=300)
plot2
dev.off()

#tiff("ScatterPlot_MitoFeatures.tiff", units="in", width=5, height=5, res=300)
#plot1 + plot2
#dev.off()


# Normalise:
#seurat.object <- NormalizeData(seurat.object)
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify the most highly variable genes:
seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 2000)

# Save these genes in tables:
top2000 <- head(VariableFeatures(seurat.object), 2000)
write.table(top2000, file="SeuratTopDEG_Top2000", quote=FALSE, sep='\t', col.names = FALSE)

top1000 <- head(VariableFeatures(seurat.object), 1000)
write.table(top1000, file="SeuratTopDEG_Top1000", quote=FALSE, sep='\t', col.names = FALSE)

top100 <- head(VariableFeatures(seurat.object), 100)
write.table(top100, file="SeuratTopDEG_Top100", quote=FALSE, sep='\t', col.names = FALSE)

# Make graphs of the top 10 and 20 genes:
top10 <- head(VariableFeatures(seurat.object), 10)
top20 <- head(VariableFeatures(seurat.object), 20)

# Plot the top (10) most highly variable genes:
plot1 <- VariableFeaturePlot(seurat.object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot3 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)

tiff("VariableGenesPlot.tiff", units="in", width=15, height=15, res=300)
plot1
dev.off()

tiff("Top10VariableGenesPlot.tiff", units="in", width=15, height=15, res=300)
plot2
dev.off()

tiff("Top20VariableGenesPlot.tiff", units="in", width=15, height=15, res=300)
plot3
dev.off()

# Scaling the data:
all.genes <- rownames(seurat.object)
seurat.object <- ScaleData(seurat.object, features = all.genes)

# Perform linear dimensional reduction:
seurat.object <- RunPCA(seurat.object, features = VariableFeatures(object = seurat.object))
print(seurat.object[["pca"]], dims = 1:10, nfeatures = 6)

# Plot PCA and Heatmaps:
tiff("DimRedPlot.tiff", units="in", width=20, height=30, res=300)
VizDimLoadings(seurat.object, dims = 1:10, reduction = "pca")
dev.off()

tiff("PCAPlot.tiff", units="in", width=7, height=5, res=300)
DimPlot(seurat.object, reduction = "pca")
dev.off()

tiff("HeatmapPlot.tiff", units="in", width=7, height=5, res=300)
par(mai=c(6,6,6,10))
DimHeatmap(seurat.object, dims = 1, cells = 200, balanced = TRUE, fast = FALSE) + ggplot2::scale_fill_gradientn(colors = c("steelblue1", "white", "tomato")) + theme_grey(base_size = 8) + theme(axis.text.x = element_blank())
dev.off()

tiff("MultiHeatmapPlot.tiff", units="in", width=40, height=30, res=300)
par(mai=c(6,6,6,10))
DimHeatmap(seurat.object, dims = 1:12, cells = 200, balanced = TRUE, fast = FALSE) + theme_grey(base_size = 8) + theme(axis.text.x = element_blank())
dev.off()


# Define dimensionality of data:
seurat.object <- JackStraw(seurat.object, num.replicate = 2000)
seurat.object <- ScoreJackStraw(seurat.object, dims = 1:20)

tiff("StrawPlot.tiff", units="in", width=5, height=5, res=300)
JackStrawPlot(seurat.object, dims = 1:20)
dev.off()

tiff("ElbowPlot.tiff", units="in", width=5, height=5, res=300)
ElbowPlot(seurat.object)
dev.off()


