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

# MAPQ scores for gmap (from: https://www.researchgate.net/post/How_to_derive_multiple_mapped_reads_from_a_SAM_file)
# 40: Maps uniquely to one location on the genome
#  3: Maps to two locations on the genome
#  2: Maps to three locations
#  1: maps to 5 - 10 locations
#  0: maps to more than 10 locations

cd $MainDirectory

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "##### Trinity Constructs Unaligned for ${id}_${project} #####" >> ${MainDirectory}/StepSummaries/SummaryOfGmapUnaligned; echo "" >> ${MainDirectory}/StepSummaries/SummaryOfGmapUnaligned; echo "## Raw Assembly:" >> ${MainDirectory}/StepSummaries/SummaryOfGmapUnaligned; grep -c "No paths found for" GmapRaw_${id}_${project}.e* >> ${MainDirectory}/StepSummaries/SummaryOfGmapUnaligned; echo "## Filtered Assembly:" >> ${MainDirectory}/StepSummaries/SummaryOfGmapUnaligned; grep -c "No paths found for" GmapFiltered_${id}_${project}.e* >> ${MainDirectory}/StepSummaries/SummaryOfGmapUnaligned; echo -e "\n********************\n" >> ${MainDirectory}/StepSummaries/SummaryOfGmapUnaligned;done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "##### MAPQ for Gmap alignment of ${id}_${project} to genome #####" >> ${MainDirectory}/StepSummaries/SummaryOfGmapMAPQ; echo "" >> ${MainDirectory}/StepSummaries/SummaryOfGmapMAPQ; echo "## Raw Assembly:" >> ${MainDirectory}/StepSummaries/SummaryOfGmapMAPQ; grep -v "^@"  ${MainDirectory}/GmapRaw_${id}_${project}/Trinity_${id}_${project}.genome.sam | cut -f5 | $StatScript >> ${MainDirectory}/StepSummaries/SummaryOfGmapMAPQ; grep -v "^@" ${MainDirectory}/GmapRaw_${id}_${project}/Trinity_${id}_${project}.genome.sam | cut -f5 | sort | uniq -c | sort -k1,1 -nbr | head -n5 >> ${MainDirectory}/StepSummaries/SummaryOfGmapMAPQ; echo "## Filtered Assembly:" >> ${MainDirectory}/StepSummaries/SummaryOfGmapMAPQ; grep -v "^@"  ${MainDirectory}/GmapFiltered_${id}_${project}/Transrate_${id}_${project}.genome.sam | cut -f5 | $StatScript >> ${MainDirectory}/StepSummaries/SummaryOfGmapMAPQ; grep -v "^@" ${MainDirectory}/GmapFiltered_${id}_${project}/Transrate_${id}_${project}.genome.sam | cut -f5 | sort | uniq -c | sort -k1,1 -nbr | head -n5 >> ${MainDirectory}/StepSummaries/SummaryOfGmapMAPQ; echo -e "\n********************\n" >> ${MainDirectory}/StepSummaries/SummaryOfGmapMAPQ;done; done

############################################
##            WRITE TO LOG FILE           ##
############################################

############################################
##          MOVE FILES TO OUTPUT          ##
############################################
rm -rf Gmap_*.pbs
mv GmapFiltered_*.pbs GmapRaw_*.pbs ${MainDirectory}/PBSin
mv Gmap*_*.e* ${MainDirectory}/PBSout
mv Gmap*_*.o* ${MainDirectory}/PBSout
cd $MainDirectory

############################################
##             STEP  FINISHED             ##
############################################


