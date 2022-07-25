
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


# WARNING!!!dims = 1:N
# --- Edit --- Replace through the file:
# Find "dims = 1:N"
# Replace by "dims = 1:X"


# Reload the object:
rsem.melange <- readRDS(file = "rsem_melange.rds")

#  --- Edit --- Create a new Seurat object using the data from melange:
# with protein-coding genes:
seurat.object <- CreateSeuratObject(counts = rsem.melange@data, project = "MySeuratProject", assay = "RNA", min.cells = 3, min.features = 900)
ncol(seurat.object)

# without protein-coding genes:
#seurat.object <- CreateSeuratObject(counts = rsem.melange@data, project = "MySeuratProject", assay = "RNA", min.cells = 3, min.features = 300)
#ncol(seurat.object)


# Identify mitochondrial reads:
seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, pattern = "^MT-")
seurat.object[["percent.mt"]][is.na(seurat.object[["percent.mt"]])] <- 0

# Normalise:
#seurat.object <- NormalizeData(seurat.object)
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify the most highly variable genes:
seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 2000)

# Scaling the data:
all.genes <- rownames(seurat.object)
seurat.object <- ScaleData(seurat.object, features = all.genes)

# Perform linear dimensional reduction:
seurat.object <- RunPCA(seurat.object, features = VariableFeatures(object = seurat.object))

# --- Edit --- Test possible resolutions:
seurat.resolution <- FindNeighbors(seurat.object, dims = 1:N)

seurat.resolution <- FindClusters(seurat.resolution, reduction.type="umap", resolution = c(0.2, 0.5, 0.8, 0.6, 1.0, 1.2), dims = 1:N, graph.name = "RNA_snn")

tiff("Clustree.tiff", units="in", width=7, height=9, res=300)
clustree(seurat.resolution)
dev.off()

# Backup Seurat object before clustering:
backup.seurat.object <- seurat.object

# --- Edit --- Cluster Cells (chosen resolution):
seurat.object <- FindNeighbors(seurat.object, dims = 1:N)
seurat.object <- FindClusters(seurat.object, resolution = 0.5)
head(Idents(seurat.object), 5)

# --- Edit --- Run clustering:
seurat.object <- RunUMAP(seurat.object, dims = 1:N)

tiff("UMAP.tiff", units="in", width=5, height=5, res=300)
# add to line below if needed +xlim(-5,10)+ylim(-6,3)
DimPlot(seurat.object, reduction = "umap")
dev.off()

# Create table of clusters:
write.table(seurat.object@active.ident, file='Cell_Clusters.tsv', quote=FALSE, sep='\t', col.names = FALSE)

# Save object as is:
saveRDS(seurat.object, file = "seurat.object_clustered.rds")

# Assign markers (to an individual cluster):
#cluster1.markers <- FindMarkers(seurat.object, ident.1 = 2, min.pct = 0.25)
#head(cluster1.markers, n = 5)

# Assign markers (to all clusters):
seurat.object.markers <- FindAllMarkers(seurat.object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.table(seurat.object.markers, file="SeuratMarkers", quote=FALSE, sep='\t', col.names = FALSE)

# --- Edit --- Cluster Cells (lower resolution):
seurat.object <- FindNeighbors(seurat.object, dims = 1:N)
seurat.object <- FindClusters(seurat.object, resolution = 0.2)
head(Idents(seurat.object), 5)

# --- Edit --- Run clustering:
seurat.object <- RunUMAP(seurat.object, dims = 1:N)

tiff("UMAP_lowerRes.tiff", units="in", width=5, height=5, res=300)
# add to line below if needed +xlim(-5,10)+ylim(-6,3)
DimPlot(seurat.object, reduction = "umap")
dev.off()

# Create table of clusters:
write.table(seurat.object@active.ident, file='Cell_Clusters_lowerRes.tsv', quote=FALSE, sep='\t', col.names = FALSE)

# --- Edit --- Cluster Cells (higher resolution):
seurat.object <- FindNeighbors(seurat.object, dims = 1:N)
seurat.object <- FindClusters(seurat.object, resolution = 0.8)
head(Idents(seurat.object), 5)

# --- Edit --- Run clustering:
seurat.object <- RunUMAP(seurat.object, dims = 1:N)

tiff("UMAP_higherRes.tiff", units="in", width=5, height=5, res=300)
# add to line below if needed +xlim(-5,10)+ylim(-6,3)
DimPlot(seurat.object, reduction = "umap")
dev.off()

# Create table of clusters:
write.table(seurat.object@active.ident, file='Cell_Clusters_higherRes.tsv', quote=FALSE, sep='\t', col.names = FALSE)

