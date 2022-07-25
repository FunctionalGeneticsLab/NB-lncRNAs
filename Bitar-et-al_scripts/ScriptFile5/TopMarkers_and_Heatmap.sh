
#!/bin/sh

# User must supply the directory name.

dir=$1; currdir=`pwd`

MarkersHeatmap=/working/lab_julietF/mainaB/Scripts/MarkersHeatmap.R

## Markers from Literature:

markerslist=/working/lab_julietF/mainaB/Project_NormalBreast/scRNAseq/LiteratureMarkers/Markers_SeuratNaming

topmarkerslist=/working/lab_julietF/mainaB/Project_NormalBreast/scRNAseq/LiteratureMarkers/Top180Markers_SeuratNaming


# Enter directory containing Seruat output files:

cd ${currdir}/${dir}/SeuratResults

nclusters=`sort -nbr -k2,2 Cell_Clusters.tsv | head -n1 | cut -f2`

echo "===> Clusters 0 to ${nclusters} were found."
echo ""

## Create lists with top 10 or 20 markers per cluster:

# Top 10:
seq 0 ${nclusters} | while read n; do cut -f7,8 SeuratMarkers | grep -Pe "^${n}\t" | head | cut -f2 >> Top10PerCluster_SeuratMarkers; done

# Top 20:
seq 0 ${nclusters} | while read n; do cut -f7,8 SeuratMarkers | grep -Pe "^${n}\t" | head -n20 | cut -f2 >> Top20PerCluster_SeuratMarkers; done


## Sort Literature Markers:

cat ${markerslist} | tr '_' '-' | while read id; do grep -Pe "\t${id}$" SeuratMarkers | cut -f7,8 >> temporary; done

sort -k1,1 -nb temporary | cut -f2 >> LiteratureMarkers_in_SeuratMarkers; rm -rf temporary

cat ${topmarkerslist} | tr '_' '-' | while read id; do grep -Pe "\t${id}$" SeuratMarkers | cut -f7,8 >> temporary; done

sort -k1,1 -nb temporary | cut -f2 >> TopLiteratureMarkers_in_SeuratMarkers; rm -rf temporary


## Run routine to make Heatmaps of genes of interest:

# Top 10:
Rscript ${MarkersHeatmap} Top10PerCluster_SeuratMarkers Top10PerCluster.tiff

# Top 20:
Rscript ${MarkersHeatmap} Top20PerCluster_SeuratMarkers Top20PerCluster.tiff

# Literature:
Rscript ${MarkersHeatmap} LiteratureMarkers_in_SeuratMarkers LiteratureMarkers.tiff

# Top 100 from Literature:
Rscript ${MarkersHeatmap} TopLiteratureMarkers_in_SeuratMarkers TopLiteratureMarkers.tiff

## Content of MarkersHeatmap.R:

#!/usr/bin/env Rscript
# module load R/4.0.2
# Load libraries:
#library(melange)
#library(dplyr)
#library(Seurat)
#library(patchwork)
#library(ggplot2)
#library(clustree)
#library(stringr)
# Load Seurat object
#seurat.object <- readRDS(file = "seurat.object_clustered.rds")
# Load gene list:
#list <- as.list(read.table(args[1]), header = FALSE)
# Create output:
#tiff(args[2], units="in", width=30, height=5, res=300)
#DoHeatmap(seurat.object, assay = "RNA", group.by= "seurat_clusters", features = list$V1) + theme(axis.text.y = element_text(size = 5))
#dev.off()

