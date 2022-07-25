
#!/usr/bin/env Rscript
# module load R/4.0.2

# Use arguments:
args<-commandArgs(TRUE)
# 1-Dir; 2-RootCluster; 3-Dimensions;

# Load libraries
library(slingshot)
library(Seurat)
library(SingleCellExperiment)
library(Polychrome)
library(ggbeeswarm)
library(ggthemes)
library(scales)
library(viridis)

# Define Function - Start
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}
# Define Function - End


### ROUTINE 0:
# Read in object created by Seurat run:
infile <- paste(args[1],"seurat.object_clustered.rds", sep="/")
seurat.object  <- readRDS(file = infile)

# Transform the object into Sligshot format:
sce <- as.SingleCellExperiment(seurat.object)
table(sce$ident)


### ROUTINE 1:
# Define lineages with automatic set of root-state:
sce1 <- slingshot(sce, clusterLabels = sce$ident, reducedDim = "UMAP", allow.breaks = FALSE)
lnes1 <- getLineages(reducedDim(sce1,"UMAP"), sce1$ident)
lnes1@lineages

# Define lineages with manual set of root-state:
sce2 <- slingshot(sce, clusterLabels = sce$ident, reducedDim = "UMAP", allow.breaks = FALSE, start.clus=args[2])
lnes2 <- getLineages(reducedDim(sce2,"UMAP"), sce2$ident, start.clus = args[2])
lnes2@lineages

# Define the cluster colors (user can change it with own color scheme):
my_color <- createPalette(length(levels(sce$ident)), c("#FF0000", "#FF8000", "#FFFF00", "#00FF00", "#0080FF", "#7F00FF", "#FF00FF"), M=1000)
names(my_color) <- unique(as.character(sce$ident))

# Create dataframe:
slingshot_df1 <- data.frame(colData(sce1))
slingshot_df2 <- data.frame(colData(sce2))

# Re-order y-axis for better figure (user would need to tailor cluster names):
#slingshot_df$ident = factor(slingshot_df$ident, levels=c(4,2,1,0,3,5,6))

# Plot clusters as Seurat defined:
tiff("Slingshot_Trajectories_NoClustering_NoRoot.tiff", units="in", width=5, height=5, res=300)
plot(reducedDims(sce1)$UMAP, col = my_color[as.character(sce1$ident)], pch=16, asp = 1)
legend("bottomleft",legend = names(my_color[levels(sce1$ident)]), fill = my_color[levels(sce1$ident)])
lines(SlingshotDataSet(lnes1), lwd=2, type = 'lineages', col = c("black"))
dev.off()

tiff("Slingshot_Trajectories_NoClustering_Root.tiff", units="in", width=5, height=5, res=300)
plot(reducedDims(sce2)$UMAP, col = my_color[as.character(sce2$ident)], pch=16, asp = 1)
legend("bottomleft",legend = names(my_color[levels(sce2$ident)]), fill = my_color[levels(sce2$ident)])
lines(SlingshotDataSet(lnes2), lwd=2, type = 'lineages', col = c("black"))
dev.off()


# Plot clusters according to pseudotime:
tiff("Slingshot_ClusterPseudotimeOrder_NoRoot.tiff", units="in", width=5, height=5, res=300)
ggplot(slingshot_df1, aes(x = slingPseudotime_1, y = ident, colour = ident)) +
    geom_quasirandom(groupOnX = FALSE) + theme_classic() +
    xlab("First Slingshot pseudotime") + ylab("cell type") +
    ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)
dev.off()

tiff("Slingshot_ClusterPseudotimeOrder_Root.tiff", units="in", width=5, height=5, res=300)
ggplot(slingshot_df2, aes(x = slingPseudotime_1, y = ident, colour = ident)) +
    geom_quasirandom(groupOnX = FALSE) + theme_classic() +
    xlab("First Slingshot pseudotime") + ylab("cell type") +
    ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)
dev.off()

# Plot clusters according to pseudotime1 x pseudotime2:
tiff("Slingshot_ClusterPseudotimeVersus_NoRoot.tiff", units="in", width=5, height=5, res=300)
ggplot(slingshot_df1, aes(x = slingPseudotime_1, y = slingPseudotime_2, colour = ident)) +
    geom_quasirandom(groupOnX = FALSE) + theme_classic() +
    xlab("First Slingshot pseudotime") + ylab("Second Slingshot pseudotime") +
    ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)
dev.off()

tiff("Slingshot_ClusterPseudotimeVersus_Root.tiff", units="in", width=5, height=5, res=300)
ggplot(slingshot_df2, aes(x = slingPseudotime_1, y = slingPseudotime_2, colour = ident)) +
    geom_quasirandom(groupOnX = FALSE) + theme_classic() +
    xlab("First Slingshot pseudotime") + ylab("Second Slingshot pseudotime") +
    ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)
dev.off()

### ROUTINE 2:
# Define cluster colours:
#cluster_colors <- cell_pal(seurat.object$seurat_clusters, brewer_pal("qual", "Set2"))
cluster_colors <- cell_pal(seurat.object$seurat_clusters, hue_pal())

# Plot Seurat clusters for control:
tiff("Slingshot_SeuratClusters.tiff", units="in", width=5, height=5, res=300)
DimPlot(seurat.object, pt.size = 0.5, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
dev.off()

seu <- seurat.object
seu <- RunUMAP(seu, dims = 1:args[3], seed.use = 4867)

seu <- FindNeighbors(seu, verbose = FALSE, dims = 1:args[3])
seu <- FindClusters(seu, algorithm = 1, random.seed = 256, resolution = 1.2)
head(Idents(seu), 5)


### ROUTINE 3A:
# Run slingshot:
sds <- slingshot(Embeddings(seu, "umap"), clusterLabels = seu$seurat_clusters, stretch = 0, start.clus=args[2])
lnes <- getLineages(reducedDim(sds,"UMAP"), sce$ident, start.clus=args[2])
lnes@lineages


# Plot trajectories:
tiff("Slingshot_Trajectory_Root.tiff", units="in", width=5, height=5, res=300)
plot(reducedDim(sds), col = cluster_colors, pch = 16, cex = 0.5)
lines(sds, lwd = 2, type = 'lineages', col = 'black')
dev.off()

tiff("Slingshot_TrajetoryCurve_Root.tiff", units="in", width=5, height=5, res=300)
plot(reducedDim(sds), col = cluster_colors, pch = 16, cex = 0.5)
lines(sds, lwd = 2, col = 'black')
dev.off()


# Plot pseudotime:
nr <- 2
nc <- 2

pt <- slingPseudotime(sds)
nms <- colnames(pt)
pal <- viridis(100, end = 0.95)
#nr <- ceiling(length(nms)/nc)

tiff("Slingshot_Pseudotime_Root.tiff", units="in", width=15, height=10, res=300)
par(mfrow = c(nr, nc))

for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(sds, lwd = 1, col = 'black', type = 'lineages')
}
dev.off()

# Print cluster labels:
write.table(sds@clusterLabels, file='CellClustersSlingshot_Root.tsv', quote=FALSE, sep='\t', col.names = FALSE)



### ROUTINE 3B:
# Run slingshot:
sdsnr <- slingshot(Embeddings(seu, "umap"), clusterLabels = seu$seurat_clusters, stretch = 0)
lnesnr <- getLineages(reducedDim(sdsnr,"UMAP"), sce$ident)
lnesnr@lineages

# Plot trajectories:
tiff("Slingshot_Trajectory_NoRoot.tiff", units="in", width=5, height=5, res=300)
plot(reducedDim(sdsnr), col = cluster_colors, pch = 16, cex = 0.5)
lines(sdsnr, lwd = 2, type = 'lineages', col = 'black')
dev.off()

tiff("Slingshot_TrajetoryCurve_NoRoot.tiff", units="in", width=5, height=5, res=300)
plot(reducedDim(sdsnr), col = cluster_colors, pch = 16, cex = 0.5)
lines(sdsnr, lwd = 2, col = 'black')
dev.off()


# Plot pseudotime:
nr <- 2
nc <- 2

pt <- slingPseudotime(sdsnr)
nms <- colnames(pt)
pal <- viridis(100, end = 0.95)
#nr <- ceiling(length(nms)/nc)

tiff("Slingshot_Pseudotime_NoRoot.tiff", units="in", width=15, height=10, res=300)
par(mfrow = c(nr, nc))

for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(sdsnr), col = colors, pch = 16, cex = 0.5, main = i)
  lines(sdsnr, lwd = 1, col = 'black', type = 'lineages')
}
dev.off()

# Print cluster labels:
write.table(sdsnr@clusterLabels, file='CellClustersSlingshot_NoRoot.tsv', quote=FALSE, sep='\t', col.names = FALSE)


# Not Used:

### BLOCK 1
# Plot clusters as Seurat defined:
#tiff("Slingshot_Lineages_NoRoot.tiff", units="in", width=5, height=5, res=300)
#plot(reducedDims(sce1)$PCA, col = my_color[as.character(sce1$ident)], pch=16, asp = 1)
#legend("bottomleft",legend = names(my_color[levels(sce1$ident)]), fill = my_color[levels(sce1$ident)])
#lines(SlingshotDataSet(lnes1), lwd=2, type = 'lineages', col = c("black"))
#dev.off()

#tiff("Slingshot_Lineages_Root.tiff", units="in", width=5, height=5, res=300)
#plot(reducedDims(sce2)$PCA, col = my_color[as.character(sce2$ident)], pch=16, asp = 1)
#legend("bottomleft",legend = names(my_color[levels(sce2$ident)]), fill = my_color[levels(sce2$ident)])
#lines(SlingshotDataSet(lnes2), lwd=2, type = 'lineages', col = c("black"))
#dev.off()


### BLOCK 2
# Run UMAP on Seurat object:
#seu <- seurat.object
#seu <- RunUMAP(seu, dims = 1:args[3], seed.use = 4867)

# Plot results for control:
#tiff("Slingshot_SeuratClusters.tiff", units="in", width=5, height=5, res=300)
#DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5, label = TRUE, repel = TRUE) + scale_color_brewer(type = "qual", palette = "Set2")
#dev.off()

# Run clustering:
#seu <- FindNeighbors(seu, verbose = FALSE, dims = 1:args[3])
#seu <- FindClusters(seu, algorithm = 1, random.seed = 256, resolution = 0.5)
#head(Idents(seu), 5)



