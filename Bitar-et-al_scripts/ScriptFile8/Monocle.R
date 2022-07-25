
# Load libraries
library(reshape)
library(Seurat)
library(monocle)
library(viridis)

# Use arguments:
args<-commandArgs(TRUE)

D <- as.integer(args[1])


# Load data from Seurat object:
Seurat.object <- readRDS(file = "seurat.object_clustered.rds")

# Separate into components:
counts.data <- as(as.matrix(Seurat.object@assays$RNA@data), 'sparseMatrix')
pheno.data <- new('AnnotatedDataFrame', data = Seurat.object@meta.data)
feature.data <- data.frame(gene_short_name = row.names(counts.data), row.names = row.names(counts.data))
feature.data <- new('AnnotatedDataFrame', data = feature.data)

# Construct CellDataSet.
HSMM <- newCellDataSet(counts.data, phenoData = pheno.data, featureData = feature.data, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())

# Load SR information from SCENT:
SRtab <- read.table("Cells_withClusters_withSR_sorted", header=T, sep="\t", row.names=1)
pData(HSMM)$SRlnc <- -1/log10(SRtab$Lnc_logSR)
pData(HSMM)$SRn13 <- -1/log10(SRtab$Net13_logSR)
pData(HSMM)$SRn17 <- -1/log10(SRtab$Net17_logSR)

# Load ZEB1 information for cells:
ZEBtab <- read.table("ZEB1_anyTPM_sorted", header=T, sep="\t", row.names=1)
pData(HSMM)$logZEB1 <- log10(ZEBtab$ZEB1+1)
pData(HSMM)$rawZEB1 <- ZEBtab$ZEB1


print(head(pData(HSMM)))

# Estimate size factors to help normalize for differences in RNA levels across cells:
HSMM <- estimateSizeFactors(HSMM)

# Estimate dispersion:
# Warning: only needed, when using a CellDataSet with a negbinomial or negbinomial.size expression 
HSMM <- estimateDispersions(HSMM)

# Remove genes with very low expression levels:
HSMM <- detectGenes(HSMM, min_expr = 0.1)

# See data structure:
print(head(fData(HSMM)))

# Prepare to look at the distribution of RNA totals across the cells:
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 100000]
#HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]

# Define expression boundaries:
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) + 2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) - 2*sd(log10(pData(HSMM)$Total_mRNAs)))

# Remove cells with either very low RNA recovery or far more RNA that the typical cell:
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound | pData(HSMM)$Total_mRNAs < upper_bound]
HSMM <- detectGenes(HSMM, min_expr = 0.1)

# Plot the distribution of RNA numbers in cells per cluster:
tiff("Monocle_RNAdistribution.tiff", units="in", width=5, height=4, res=300)
qplot(Total_mRNAs, data = pData(HSMM), color = seurat_clusters, geom = "density") + geom_vline(xintercept = lower_bound) + geom_vline(xintercept = upper_bound) 
dev.off()

# Define a subset of expressed genes:
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))

# Standardize each gene, so that they are all on the same scale and melt the data with plyr for easier plotting:
L <- log(exprs(HSMM[expressed_genes,]))
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
tiff("Monocle_RNAexpression.tiff", units="in", width=5, height=4, res=300)
qplot(value, geom = "density", data = melted_dens_df) + stat_function(fun = dnorm, size = 0.5, color = 'red') + xlab("Standardized log(FPKM)") + ylab("Density")
dev.off()
    
# Run dimensionality reduction:
# max_components = the number of dimensions to be returned.
# num_dim = the number of dimension for the PCA preprocessing.
HSMM <- reduceDimension(HSMM, max_components = D, method = 'DDRTree')

# Order cells along the inferred trajectory:
HSMM <- orderCells(HSMM)

# Plot trajectory:
tiff("Monocle_Trajectory_ColByCluster.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")
dev.off()

tiff("Monocle_Trajectory_Pseudotime.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
dev.off()

# Plot branches separately:
tiff("Monocle_Trajectory_Branches.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "seurat_clusters") + facet_wrap(~State, nrow = 1)
dev.off()

# Plot ZEB1 expression:
tiff("Monocle_Trajectory_ZEB1_ColourC.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "logZEB1") + scale_colour_viridis(limits = c(min(HSMM$logZEB1, na.rm = TRUE),max(HSMM$logZEB1, na.rm = TRUE)), direction = 1, option = "C", na.value = "yellow"); dev.off()
tiff("Monocle_Trajectory_ZEB1_ColourD.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "logZEB1") + scale_colour_viridis(limits = c(min(HSMM$logZEB1, na.rm = TRUE),max(HSMM$logZEB1, na.rm = TRUE)), direction = 1, option = "D", na.value = "yellow"); dev.off()


# Plot cells colured by SR values:
b <- seq.int(min(HSMM$SRlnc, na.rm = TRUE),max(HSMM$SRlnc, na.rm = TRUE),0.5)
b <- format(round(b, 2), nsmall=2)
b <- as.numeric(b)
#m=median(HSMM$SRlnc)
m <- min(HSMM$SRlnc, na.rm = TRUE)+((max(HSMM$SRlnc, na.rm = TRUE)-min(HSMM$SRlnc, na.rm = TRUE))/2)
tiff("Monocle_Trajectory_SRLnc_ColourA.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "SRlnc") + scale_color_gradient2(breaks=b, low="#6aa84f", high="#674ea7", guide = "colourbar", midpoint=m); dev.off()
tiff("Monocle_Trajectory_SRLnc_ColourB.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "SRlnc") + scale_color_gradient2(breaks=b, low="#ead1dc", mid="#d5a6bd", high="#c90076", guide = "colourbar", midpoint=m); dev.off()
tiff("Monocle_Trajectory_SRLnc_ColourC.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "SRlnc") + scale_colour_viridis(limits = c(min(HSMM$SRlnc, na.rm = TRUE),max(HSMM$SRlnc, na.rm = TRUE)-1), direction = 1, option = "C", na.value = "yellow"); dev.off()
tiff("Monocle_Trajectory_SRLnc_ColourD.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "SRlnc") + scale_colour_viridis(limits = c(min(HSMM$SRlnc, na.rm = TRUE),max(HSMM$SRlnc, na.rm = TRUE)-1), direction = 1, option = "D", na.value = "yellow"); dev.off()


b <- seq.int(min(HSMM$SRn17, na.rm = TRUE),max(HSMM$SRn17, na.rm = TRUE),0.5)
b <- format(round(b, 2), nsmall=2)
b <- as.numeric(b)
#m=median(HSMM$SRn17)
m <- min(HSMM$SRn17, na.rm = TRUE)+((max(HSMM$SRn17, na.rm = TRUE)-min(HSMM$SRn17, na.rm = TRUE))/2)
tiff("Monocle_Trajectory_SRNet17_ColourA.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "SRn17") + scale_color_gradient2(breaks=b, low="#6aa84f", high="#674ea7", guide = "colourbar", midpoint=m); dev.off()
tiff("Monocle_Trajectory_SRNet17_ColourB.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "SRn17") + scale_color_gradient2(breaks=b, low="#ead1dc", mid="#d5a6bd", high="#c90076", guide = "colourbar", midpoint=m); dev.off()

tiff("Monocle_Trajectory_SRNet17_ColourC.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "SRn17") + scale_colour_viridis(limits = c(min(HSMM$SRn17, na.rm = TRUE),max(HSMM$SRn17, na.rm = TRUE)-1), direction = 1, option = "C", na.value = "yellow"); dev.off()
tiff("Monocle_Trajectory_SRNet17_ColourD.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "SRn17") + scale_colour_viridis(limits = c(min(HSMM$SRn17, na.rm = TRUE),max(HSMM$SRn17, na.rm = TRUE)-1), direction = 1, option = "D", na.value = "yellow"); dev.off()


# Function for root state:
root_lnc <- as.numeric(HSMM$State[which.max(HSMM$SRlnc)])
root_n13 <- as.numeric(HSMM$State[which.max(HSMM$SRn13)])
root_n17 <- as.numeric(HSMM$State[which.max(HSMM$SRn17)])
root_lnc; root_n13; root_n17
#root <- root_lnc

# Order cells along the inferred trajectory, defining root state:
HSMM <- orderCells(HSMM, root_state = root_lnc)

# Plot rooted trajectory:
tiff("Monocle_Trajectory_ColByCluster_RootLnc.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "seurat_clusters"); dev.off()

tiff("Monocle_Trajectory_Pseudotime_RootLnc.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "Pseudotime"); dev.off()

# Plot rooted trajectory coloured by SR values:
b <- seq.int(min(HSMM$SRlnc, na.rm = TRUE),max(HSMM$SRlnc, na.rm = TRUE),0.5)
b <- format(round(b, 2), nsmall=2)
b <- as.numeric(b)
#m=median(HSMM$SRlnc)
m <- min(HSMM$SRlnc, na.rm = TRUE)+((max(HSMM$SRlnc, na.rm = TRUE)-min(HSMM$SRlnc, na.rm = TRUE))/2)
tiff("Monocle_Trajectory_SRLnc_RootLnc_ColourA.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "SRlnc") + scale_color_gradient2(breaks=b, low="#6aa84f", high="#674ea7", guide = "colourbar", midpoint=m); dev.off()
tiff("Monocle_Trajectory_SRLnc_RootLnc_ColourB.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "SRlnc") + scale_color_gradient2(breaks=b, low="#ead1dc", mid="#d5a6bd", high="#c90076", guide = "colourbar", midpoint=m); dev.off()
tiff("Monocle_Trajectory_SRLnc_RootLnc_ColourC.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "SRlnc") + scale_colour_viridis(limits = c(min(HSMM$SRlnc, na.rm = TRUE),max(HSMM$SRlnc, na.rm = TRUE)-1), direction = 1, option = "C", na.value = "yellow"); dev.off()
tiff("Monocle_Trajectory_SRLnc_ColourD.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "SRlnc") + scale_colour_viridis(limits = c(min(HSMM$SRlnc, na.rm = TRUE),max(HSMM$SRlnc, na.rm = TRUE)-1), direction = 1, option = "D", na.value = "yellow"); dev.off()


# Order cells along the inferred trajectory, defining root state:
HSMM <- orderCells(HSMM, root_state = root_n17)

# Plot rooted trajectory:
tiff("Monocle_Trajectory_ColByCluster_RootN17.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "seurat_clusters"); dev.off()

tiff("Monocle_Trajectory_Pseudotime_RootN17.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "Pseudotime"); dev.off()

# Plot rooted trajectory coloured by SR values:
b <- seq.int(min(HSMM$SRn17, na.rm = TRUE),max(HSMM$SRn17, na.rm = TRUE),0.5)
b <- format(round(b, 2), nsmall=2)
b <- as.numeric(b)
#m=median(HSMM$SRn17)
m <- min(HSMM$SRn17, na.rm = TRUE)+((max(HSMM$SRn17, na.rm = TRUE)-min(HSMM$SRn17, na.rm = TRUE))/2)
tiff("Monocle_Trajectory_SRnet17_RootNet17_ColourA.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "SRn17") + scale_color_gradient2(breaks=b, low="#6aa84f", high="#674ea7", guide = "colourbar", midpoint=m); dev.off()
tiff("Monocle_Trajectory_SRnet17_RootNet17_ColourB.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "SRn17") + scale_color_gradient2(breaks=b, low="#ead1dc", mid="#d5a6bd", high="#c90076", guide = "colourbar", midpoint=m); dev.off()
tiff("Monocle_Trajectory_SRnet17_RootNet17_ColourC.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "SRn17") + scale_colour_viridis(limits = c(min(HSMM$SRn17, na.rm = TRUE),max(HSMM$SRn17, na.rm = TRUE)-1), direction = 1, option = "C", na.value = "yellow"); dev.off()
tiff("Monocle_Trajectory_SRnet17_ColourD.tiff", units="in", width=5, height=4, res=300)
plot_cell_trajectory(HSMM, color_by = "SRn17") + scale_colour_viridis(limits = c(min(HSMM$SRn17, na.rm = TRUE),max(HSMM$SRn17, na.rm = TRUE)-1), direction = 1, option = "D", na.value = "yellow"); dev.off()


