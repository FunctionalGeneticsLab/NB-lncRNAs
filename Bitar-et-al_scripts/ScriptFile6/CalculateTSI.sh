#!/bin/sh

# User must supply the directory name.

dir=$1; currdir=`pwd`

markerslist=/working/lab_julietF/mainaB/scRNAseq/LiteratureMarkers/Markers_Revised_Annotation

AverageRows=/working/lab_julietF/mainaB/Scripts/AverageRows.R

MakeMDS=/working/lab_julietF/mainaB/Scripts/MakeMDS.R


# Enter input directory and create output directory

cd ${currdir}/${dir}

mkdir ${currdir}/${dir}/TSPEXinputs


### Formatting the count matrix.

# A count matrix as an input for TSPEX should have the clusters in column and features in row.

# For each cell cluster, combine all gene-level RSEM count files into a matrix:

cp SeuratResults/Cell_Clusters.tsv ${currdir}/${dir}

nclusters=`sort -nbr -k2,2 Cell_Clusters.tsv | head -n1 | cut -f2`

echo "===> Clusters 0 to ${nclusters} were found."
echo ""

### LOOP to create the integrated matrix per cluster:

seq 0 ${nclusters} | while read n; do echo "==> Entering the main step for cluster ${n}..."; clustersize=`grep -c "${n}$" Cell_Clusters.tsv`; nonseed=`echo "${clustersize}-2" | bc`; grep "${n}$" Cell_Clusters.tsv | cut -f1 | head -n2 | head -n1 | while read cellid; do cut -f1,6 BowtieRsemExpressionMatrices/RSEM_${cellid}.genes.results | sed "s/\tTPM/\t${cellid}/g" >> temporaryA; done; grep "${n}$" Cell_Clusters.tsv | cut -f1 | head -n2 | tail -n1 | while read cellid; do cut -f1,6 BowtieRsemExpressionMatrices/RSEM_${cellid}.genes.results | sed "s/\tTPM/\t${cellid}/g" >> temporaryB; done ; join -t$'\t' temporaryA temporaryB >> SeedJoin.genes.results; rm -rf temporary*; grep "${n}$" Cell_Clusters.tsv | cut -f1 | tail -n${nonseed} | while read cellid; do cut -f1,6 BowtieRsemExpressionMatrices/RSEM_${cellid}.genes.results | sed "s/\tTPM/\t${cellid}/g" >> temporary; join -t$'\t' temporary SeedJoin.genes.results >> CurrentJoin.genes.results; rm -rf SeedJoin.genes.results; mv CurrentJoin.genes.results SeedJoin.genes.results; echo "ok" >> finished; rm -rf temporary; done; rm -rf finished; count=`head -n1 SeedJoin.genes.results | tr '\t' '\n' | wc -l`; ncount=`echo "${count}-1" | bc`; echo "Cluster ${n} has ${ncount} cells."; mv SeedJoin.genes.results TSPEXinputs/RSEMmatrix_Cluster${n}; echo "Finished the main step for cluster ${n}!"; done

echo "Done!"; echo ""

### Calculate the average expression per gene for each cluster:

# Enter directory:
cd TSPEXinputs

echo "===> Calculating average per cluster with ${AverageRows}..."


# Run AverageRows.R for each cluster:
module load R/3.6.2
seq 0 ${nclusters} | while read n; do Rscript ${AverageRows} RSEMmatrix_Cluster${n} AverageCounts_Cluster${n}; done

echo "Done!"; echo ""

echo "===> Consolidating Tables..."


### Consolidate the table of counts:

join -t$'\t' RSEMmatrix_Cluster1 RSEMmatrix_Cluster0 >> SeedJoin.genes.results

seq 2 ${nclusters} | while read n; do cat RSEMmatrix_Cluster${n} >> temporary; join -t$'\t' temporary SeedJoin.genes.results >> CurrentJoin.genes.results; rm -rf SeedJoin.genes.results; mv CurrentJoin.genes.results SeedJoin.genes.results; rm -rf temporary; done; mv SeedJoin.genes.results RSEM_AllClustersCount.matrix

cat RSEM_AllClustersCount.matrix | tr '\t' ',' >> RSEM_AllClustersCount.csv


### Consolidate the table of averages:

join -t$'\t' AverageCounts_Cluster1 AverageCounts_Cluster0 | sed "s/\trow_mean\trow_mean/\tCluster1\tCluster0/g" >> SeedJoin.genes.results

seq 2 ${nclusters} | while read n; do cat AverageCounts_Cluster${n} | sed "s/\trow_mean/\tCluster${n}/g" >> temporary; join -t$'\t' temporary SeedJoin.genes.results >> CurrentJoin.genes.results; rm -rf SeedJoin.genes.results; mv CurrentJoin.genes.results SeedJoin.genes.results; rm -rf temporary; done; mv SeedJoin.genes.results AllClustersAverage.matrix

cat AllClustersAverage.matrix | tr '\t' ',' >> AllClustersAverage.csv

echo "Done!"; echo ""


### Retrieve literature marker genes:

echo "Using markers defined in file: ${markerslist}."

head -n1 AllClustersAverage.csv >> Markers_AllClustersAverage.csv

cut -f1 ${markerslist} | while read id; do grep "${id}." AllClustersAverage.csv >> Markers_AllClustersAverage.csv ; done


### Retrieve Seurat marker genes:

markerseurat=${currdir}/${dir}/SeuratResults/SeuratMarkers

echo "Using markers defined in file: ${markerseurat}."

head -n1 AllClustersAverage.csv >> SeuratMarkers_AllClustersAverage.csv

cut -f8 ${markerseurat} | sort | uniq | tr '-' '_' | while read id; do grep "${id}," AllClustersAverage.csv >> SeuratMarkers_AllClustersAverage.csv ; done



### Run tspex

echo "===> Running tspex..."

module load python/3.6.1

tspex AllClustersAverage.csv TSI_AllGenes.tsv tsi

tspex Markers_AllClustersAverage.csv TSI_MarkerGenes.tsv tsi

tspex SeuratMarkers_AllClustersAverage.csv TSI_SeuratMarkerGenes.tsv tsi

echo "Done!"; echo ""

### Make MDS plots:

echo "===> Making MDS plots..."

module load R/3.6.2

groups=`echo "${nclusters}+1" | bc`

Rscript ${MakeMDS} ${groups}

echo "Done!"; echo ""

#rm -rf ${currdir}/${dir}/Cell_Clusters.tsv

echo "FINISHED"

### Content of AverageRows.R
#!/usr/bin/env Rscript
# module load R/3.6.2
#args<-commandArgs(TRUE)
#data <- read.table(args[1], sep="\t", header=TRUE, row.names="gene_id")
#data$row_mean <- rowMeans(data)
#newdata <- cbind(gene_id=rownames(data),data)
#write.table(newdata[,c("gene_id","row_mean")], file=args[2], quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


### Content of Make MDS:
#!/usr/bin/env Rscript
# module load R/3.6.2
# library(edgeR)
# args<-commandArgs(TRUE)
#
# x <- read.table("TSI_AllGenes.tsv", header = TRUE, row.names = 1, sep = "\t")
# group_v <- seq(1, args[1])
#group <- factor(group_v)
#y <- DGEList
#y <- DGEList(counts=x,group=group)
#tiff(filename = "MDSofTSI_AllGenes.tiff")
#plotMDS(y, labels=colnames(x))
#dev.off()
#
#x <- read.table("TSI_MarkerGenes.tsv", header = TRUE, row.names = 1, sep = "\t")
#group_v <- seq(1, args[1])
#group <- factor(group_v)
#y <- DGEList
#y <- DGEList(counts=x,group=group)
#tiff(filename = "MDSofTSI_MarkerGenes.tiff")
#plotMDS(y, labels=colnames(x))
#dev.off()
#
#x <- read.table("TSI_SeuratMarkerGenes.tsv", header = TRUE, row.names = 1, sep = "\t")
#group_v <- seq(1, args[1])
#group <- factor(group_v)
#y <- DGEList
#y <- DGEList(counts=x,group=group)
#
#tiff(filename = "MDSofTSI_SeuratMarkerGenes.tiff")
#plotMDS(y, labels=colnames(x))
#dev.off()

