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

ls ${MainDirectory}/Trinity/TrinityStats_* | while read file; do sep=$(printf "%-50s" "*"); echo "${sep// /*}" >> ${MainDirectory}/StepSummaries/SummaryOfTrinityStatistics; echo "${file}" | rev | cut -d'_' -f1-2 | rev >> ${MainDirectory}/StepSummaries/SummaryOfTrinityStatistics ; grep -A1 '###\|N50\|Total\|Average' ${file} | grep -v "###\|^$\|--" | sed 's/^Stats/## Stats/g' >> ${MainDirectory}/StepSummaries/SummaryOfTrinityStatistics; done


############################################
##          MOVE FILES TO OUTPUT          ##
############################################
cd $MainDirectory
mv *Stats3nity*.pbs $MainDirectory/PBSin
mv *Stats3nity*.e* $MainDirectory/PBSout
mv *Stats3nity*.o* $MainDirectory/PBSout


############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""

