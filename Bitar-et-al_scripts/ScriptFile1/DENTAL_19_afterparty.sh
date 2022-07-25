#!/bin/sh


######################################################
#                                                    #
# ** DENTAL (de Novo Transcriptome Assembly Line) ** #
#                                                    #
######################################################
######################################################
# #             Version 2.0 - Mar 2022             # #
######################################################

############################################
##         CHECK STEP CONCLUSION!         ##
############################################

############################################
##         SET MAIN WORK DIRECTORY        ##
############################################
# Define directory:

if [ "$#" -gt 0 ]; then MainDirectory=$1; fi

if [ "$#" -eq 0 ]; then MainDirectory=`pwd`; fi

echo ""; echo "WARNING: The main directory for this run was set to ${MainDirectory}."; echo ""

############################################
##              MAKE SUMMARY              ##
############################################

cd $MainDirectory

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "##### CDhit Clustering results for ${id}_${project} #####" >> ${MainDirectory}/StepSummaries/SummaryOfCDhit; echo "# -c indicates the identity cutoff used in each run (if multiple)" >> ${MainDirectory}/StepSummaries/SummaryOfCDhit; grep -v "^\..." Filter_${id}_${project}.o* | grep " -c\|total seq\|longest and shortest\|finished" >> ${MainDirectory}/StepSummaries/SummaryOfCDhit; done; done

############################################
##          MOVE FILES TO OUTPUT          ##
############################################
mv Filter_*.pbs ${MainDirectory}/PBSin
mv *Filter_*.e* ${MainDirectory}/PBSout
mv *Filter_*.o* ${MainDirectory}/PBSout
cd $MainDirectory

############################################
##             STEP  FINISHED             ##
############################################



