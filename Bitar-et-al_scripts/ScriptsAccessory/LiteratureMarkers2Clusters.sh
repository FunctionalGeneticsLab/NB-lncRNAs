#!/bin/sh

dir=$1; currdir=`pwd`

cd ${currdir}/${dir}

n=`sort -nbr -k2,2 Cell_Clusters.tsv | head -n1 | cut -f2`

echo "===> Clusters 0 to ${n} were found."

markers="/working/lab_julietF/mainaB/Project_NormalBreast/scRNAseq/LiteratureMarkers/Markers_Revised_Annotation"

echo "Using markers in file: /working/lab_julietF/mainaB/Project_NormalBreast/scRNAseq/LiteratureMarkers/Markers_Revised_Annotation"
echo "Example of lines in the file:"
head -n3 ${markers}
echo

seq 0 ${n} | while read cluster; do echo -e "ClusterN\tMarker" >> LiteratureMarkersInCluster${cluster}; cut -f1 ${markers} | sort | uniq | while read geneid; do grep "${geneid}." SeuratMarkers | cut -f7,8 | grep -Pe "^${cluster}\t" >> LiteratureMarkersInCluster${cluster}; done; echo ${cluster}; done

echo Finished!

