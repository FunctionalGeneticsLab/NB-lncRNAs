
#!/bin/sh

# User must supply the directory name.

dir=$1; currdir=`pwd`

ref=$2

# Retrieve number of clusters:

nclusters=`sort -nbr -k2,2 ${currdir}/${dir}/SeuratResults/Cell_Clusters.tsv | head -n1 | cut -f2`

echo "===> Clusters 0 to ${nclusters} were found."
echo ""


# Enter input directory and create output directory

cd ${currdir}/${dir}

mkdir -p ${currdir}/${dir}/CompareClusters


# Retrieve reference clusters:

nrefclusters=`sort -nbr -k2,2 ${currdir}/${ref}/SeuratResults/Cell_Clusters.tsv | head -n1 | cut -f2`

seq 0 ${nrefclusters} | while read n; do grep "${n}$" ${currdir}/${ref}/SeuratResults/Cell_Clusters.tsv | cut -f1 | while read id; do grep "^${id}" ${currdir}/${dir}/SeuratResults/Cell_Clusters.tsv >> ${currdir}/${dir}/CompareClusters/CellsFrom_ReferenceCluster${n}; done; done

# Make a report:

num=`echo "${nrefclusters}+1" | bc`
echo "Reference was set to: ${ref}, which has ${num} clusters." >> ${currdir}/${dir}/CompareClusters/SummaryOfClusterDistribution
echo "This is how cells from those clusters in the reference are distributed in this set:" >> ${currdir}/${dir}/CompareClusters/SummaryOfClusterDistribution

seq 0 ${nrefclusters} | while read n; do echo "RefCluster ${n}:" >> ${currdir}/${dir}/CompareClusters/SummaryOfClusterDistribution; cut -f2 ${currdir}/${dir}/CompareClusters/CellsFrom_ReferenceCluster${n} | sort | uniq -c | sort -k1,1 -nbr | sed -e 's/^[ \t]*//' | tr ' ' '\t' >> ${currdir}/${dir}/CompareClusters/SummaryOfClusterDistribution; done


echo "Done!"; echo ""


