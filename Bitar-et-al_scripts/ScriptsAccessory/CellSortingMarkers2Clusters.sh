#!/bin/sh

dir=$1; currdir=`pwd`

cd ${currdir}/${dir}

echo "Running marker assignment for folder ${dir}..."

parent="/working/lab_julietF/mainaB/scRNAseq"

cut -f1 Cell_Clusters.tsv | sort | uniq | while read cell; do grep "$cell" ${parent}/SRAinformation/SraRunInfo_PRJNA450* | cut -d',' -f11 | while read id; do grep "$id" ${parent}/SRAinformation/SraSummary_PRJNA450* | cut -d';' -f1 >> CellTypes_raw; done; done

cut -d'"' -f2 CellTypes_raw >> CellTypes_raw1A
cut -d' ' -f2 CellTypes_raw >> CellTypes_raw1B
paste CellTypes_raw1A CellTypes_raw1B >> CellTypes_raw1

cut -f1 Cell_Clusters.tsv | sort | uniq | while read cell; do grep "$cell" ${parent}/SRAinformation/SraRunInfo_PRJNA450* | cut -d',' -f11 >> CellTypes_raw2A; done
cut -f1 Cell_Clusters.tsv | sort | uniq | while read cell; do grep "$cell" ${parent}/SRAinformation/SraRunInfo_PRJNA450* | cut -d',' -f1 >> CellTypes_raw2B; done
paste CellTypes_raw2A CellTypes_raw2B >> CellTypes_raw2

join -t$'\t' CellTypes_raw1 CellTypes_raw2 | sed "s;${parent}/SRAinformation/SraRunInfo_;;g" | sed 's/.csv:/\t/g' >> CellTypes_SeuratPass

rm -rf CellTypes_raw*

grep "_BAS_" CellTypes_SeuratPass | cut -f4 | while read cell; do grep "$cell" Cell_Clusters.tsv >> Basal_Cell_Clusters.tsv ; done

grep "_LUM_" CellTypes_SeuratPass | cut -f4 | while read cell; do grep "$cell" Cell_Clusters.tsv >> Luminal_Cell_Clusters.tsv ; done

echo "Number of Basal cells in each cluster:"

cut -f2 Basal_Cell_Clusters.tsv | sort | uniq -c

echo ""
echo "Number of Luminal cells in each cluster:"

cut -f2 Luminal_Cell_Clusters.tsv | sort | uniq -c

echo "Finished!"
echo ""

