
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
# Use arguments:
args<-commandArgs(TRUE)
# Load Seurat object
seurat.object <- readRDS(file = "seurat.object_clustered.rds")
# Load gene list:
list <- as.list(read.table(args[1]), header = FALSE)
# Create output:
tiff(args[2], units="in", width=30, height=5, res=300)
DoHeatmap(seurat.object, assay = "RNA", group.by= "seurat_clusters", features = list$V1) + theme(axis.text.y = element_text(size = 5))
dev.off()

