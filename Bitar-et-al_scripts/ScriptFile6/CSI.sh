
#!/bin/sh

# User must supply the directory name.

dir=$1; currdir=`pwd`

top=$2

# Retrieve number of clusters:

nclusters=`sort -nbr -k2,2 ${currdir}/${dir}/SeuratResults/Cell_Clusters.tsv | head -n1 | cut -f2`

echo "===> Clusters 0 to ${nclusters} were found."
echo ""


# Enter input directory and create output directory

cd ${currdir}/${dir}

mkdir -p ${currdir}/${dir}/CSI


### Retrieve the top markers.

# A count matrix as an input for TSPEX should have the clusters in column and features in row.

echo "User has defined the top "$top" markers should be used per cluster."; echo ""
echo "Generating the count matrix..."; echo ""

cd ${currdir}/${dir}/SeuratResults

#seq 0 ${nclusters} | while read n; do cut -f7,8 SeuratMarkers | grep -Pe "^${n}\t" | head -n${top} | cut -f2 >> Top${top}PerCluster_SeuratMarkers; done

seq 0 ${nclusters} | while read n; do grep -Pe "\t${n}\tE|\t${n}\tT" SeuratMarkers | awk 'NR>1{$9=$4-$5}{print}' | tr ' ' '\t' | sort -k9,9 -nbr | head -n${top} >> ${currdir}/${dir}/CSI/Top${top}ForCluster${n}_SeuratMarkers.table; done

seq 0 ${nclusters} | while read n; do cut -f8 ${currdir}/${dir}/CSI/Top${top}ForCluster${n}_SeuratMarkers.table >> ${currdir}/${dir}/CSI/Top${top}ForCluster${n}_SeuratMarkers; done

### Generate the count matrix.

cd ${currdir}/${dir}/CSI

seq 0 ${nclusters} | while read n; do head -n1 ${currdir}/${dir}/TSPEXinputs/AllClustersAverage.matrix >> Top${top}MarkersForCluster${n}_AllClustersAverage.matrix; cp Top${top}ForCluster${n}_SeuratMarkers temporary; grep "ENSG" temporary | cut -d'-' -f1 >> tempcol1; grep "ENSG" temporary | cut -d'-' -f2- >> tempcol2; grep "TRINITY" temporary | cut -d'-' -f1 >> tempcol1; grep "TRINITY" temporary | cut -d'-' -f2- | tr '-' '_' >> tempcol2; paste -d'_' tempcol1 tempcol2 | while read id; do grep -Pe "${id}\t" ${currdir}/${dir}/TSPEXinputs/AllClustersAverage.matrix >> Top${top}MarkersForCluster${n}_AllClustersAverage.matrix; done; rm -rf temporary tempcol*; done

### Run tspex

cd ${currdir}/${dir}/CSI

module load python/3.6.1

echo "===> Running tspex for Shannon method..."; seq 0 ${nclusters} | while read n; do input="${currdir}/${dir}/CSI/Top${top}MarkersForCluster${n}_AllClustersAverage.matrix"; output="CSI_Top${top}MarkersForCluster${n}"; tspex $input Shannon$output shannon_specificity; done

echo "===> Running tspex for Gini method..."; seq 0 ${nclusters} | while read n; do input="${currdir}/${dir}/CSI/Top${top}MarkersForCluster${n}_AllClustersAverage.matrix"; output="CSI_Top${top}MarkersForCluster${n}"; tspex $input Gini$output gini; done

echo "===> Running tspex for Roku method..."; seq 0 ${nclusters} | while read n; do input="${currdir}/${dir}/CSI/Top${top}MarkersForCluster${n}_AllClustersAverage.matrix"; output="CSI_Top${top}MarkersForCluster${n}"; tspex $input Roku$output roku_specificity; done

echo "===> Running tspex for Simpson method..."; seq 0 ${nclusters} | while read n; do input="${currdir}/${dir}/CSI/Top${top}MarkersForCluster${n}_AllClustersAverage.matrix"; output="CSI_Top${top}MarkersForCluster${n}"; tspex $input Simpson$output simpson; done


### Consolidate:

echo "Consolidating the TSPEX results in tables:"

stats=/working/lab_julietF/mainaB/Scripts/Statistics.sh

echo -e "Cluster\tSum\tAverage\tMedian\tMin\tMax" >> ShannonCSI_Top${top}Markers.table; seq 0 ${nclusters} | while read n; do statistics=`cut -f2 ShannonCSI_Top${top}MarkersForCluster${n} | ${stats} | cut -f2-`; echo -e "${n}\t${statistics}" >> ShannonCSI_Top${top}Markers.table; done 

echo -e "Cluster\tSum\tAverage\tMedian\tMin\tMax" >> GiniCSI_Top${top}Markers.table; seq 0 ${nclusters} | while read n; do statistics=`cut -f2 GiniCSI_Top${top}MarkersForCluster${n} | ${stats} | cut -f2-`; echo -e "${n}\t${statistics}" >> GiniCSI_Top${top}Markers.table; done 

echo -e "Cluster\tSum\tAverage\tMedian\tMin\tMax" >> RokuCSI_Top${top}Markers.table; seq 0 ${nclusters} | while read n; do statistics=`cut -f2 RokuCSI_Top${top}MarkersForCluster${n} | ${stats} | cut -f2-`; echo -e "${n}\t${statistics}" >> RokuCSI_Top${top}Markers.table; done 

echo -e "Cluster\tSum\tAverage\tMedian\tMin\tMax" >> SimpsonCSI_Top${top}Markers.table; seq 0 ${nclusters} | while read n; do statistics=`cut -f2 SimpsonCSI_Top${top}MarkersForCluster${n} | ${stats} | cut -f2-`; echo -e "${n}\t${statistics}" >> SimpsonCSI_Top${top}Markers.table; done 


echo "Done!"; echo ""
