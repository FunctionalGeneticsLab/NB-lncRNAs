#!/bin/sh

# User must supply the directory name.

dir=$2; currdir=`pwd`

module load R/4.0.2

mkdir -p ${currdir}/${dir}/GOIplots

plotoneS=/working/lab_julietF/mainaB/Scripts/PlotOneMarker_Saturated.R
plottwoS=/working/lab_julietF/mainaB/Scripts/PlotTwoMarkers_Saturated.R

plotone=/working/lab_julietF/mainaB/Scripts/PlotOneMarker.R
plottwo=/working/lab_julietF/mainaB/Scripts/PlotTwoMarkers.R

## Enter the correct directory:

cd ${currdir}/${dir}/SeuratResults

## Run if first argument is S(aturated):

if [ "$#" -lt 4 ] && [ "$1" == "S" ]; then
        plot=${plotoneS};
	Rscript ${plot} $3 ${3}.tiff
fi

if [ "$#" -gt 3 ] && [ "$1" == "S" ]; then
	plot=${plottwoS};
	Rscript ${plot} $3 $4 ${3}_and_${4}.tiff
fi

## Run if first argument is not S(aturated):

if [ "$#" -lt 4 ] && [ "$1" != "S" ]; then
        plot=${plotone};
        Rscript ${plot} $3 ${3}.tiff
fi

if [ "$#" -gt 3 ] && [ "$1" != "S" ]; then
        plot=${plottwo};
        Rscript ${plot} $3 $4 ${3}_and_${4}.tiff
fi



## Move files:

mv ${3}*.tiff ${currdir}/${dir}/GOIplots 
