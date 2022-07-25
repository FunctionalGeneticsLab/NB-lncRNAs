

#!/bin/sh

# User must supply the directory name.

dir=$1; currdir=`pwd`

ref=$2

stats="/working/lab_julietF/mainaB/Scripts/Statistics.sh"
ttestr="/working/lab_julietF/mainaB/Scripts/Ttest.R"


# Retrieve number of clusters:

nclusters=`sort -nbr -k2,2 ${currdir}/${dir}/SeuratResults/Cell_Clusters.tsv | head -n1 | cut -f2`

echo "===> Clusters 0 to ${nclusters} were found."
echo ""

# Enter input directory and create output directory

cd ${currdir}/${dir}

mkdir -p ${currdir}/${dir}/ClusterQuality

# Create total TSI table:
lastcol=`echo "${nclusters}+2" | bc`

echo -e "Cluster\tNumberOfGenes\tSum\tAverage\tMedian\tMinimum\tMaximum" >> ${currdir}/${dir}/ClusterQuality/ClusterTSI_AllGenes; i=0; seq 2 ${lastcol} | while read n; do cl=`echo "${nclusters}-${i}" | bc`; statistics=`cut -f${n} ${currdir}/${dir}/TSPEXinputs/TSI_AllGenes.tsv | grep -v Cluster | ${stats}`; echo -e "${cl}\t${statistics}" >> ${currdir}/${dir}/ClusterQuality/ClusterTSI_AllGenes; i=$((i+1)); done

# Get cluster correspondence:
cat ${currdir}/${dir}/CompareClusters/SummaryOfClusterDistribution | grep -A1 "RefC" | tr ':' '\t' | tr ' ' '\t' | tr '\n' '\t' | tr 'R' '\n' | grep efC | sed "s/^\t//g" | cut -f2,5 >> ${currdir}/${dir}/ClusterQuality/ClusterCorrespondence4ttest

cut -f1 ${currdir}/${dir}/ClusterQuality/ClusterCorrespondence4ttest | while read cnum; do grep "^${cnum}" ${currdir}/${ref}/ClusterQuality/ClusterTSI_AllGenes | cut -f4 >> tempC1; done;
cut -f2 ${currdir}/${dir}/ClusterQuality/ClusterCorrespondence4ttest | while read cnum; do grep "^${cnum}" ${currdir}/${dir}/ClusterQuality/ClusterTSI_AllGenes | cut -f4 >> tempC2; done
paste tempC1 tempC2 >> ${currdir}/${dir}/ClusterQuality/Input4ttest; rm -rf tempC*


# Run t-test routine:
cd ${currdir}/${dir}/ClusterQuality
Rscript ${ttestr} Input4ttest $3 $4 >> TtestOutput


