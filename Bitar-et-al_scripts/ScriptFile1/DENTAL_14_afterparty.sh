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
# Summarize Read Support results:
cd ${MainDirectory}

echo "** Bowtie Results **" >> ${MainDirectory}/Pipeline_Log

ls *Bowtie3nity*.e* | while read file; do echo "" >> ${MainDirectory}/Pipeline_Log; echo "$file" | cut -d'.' -f1 >> ${MainDirectory}/Pipeline_Log ; cat "$file" >> ${MainDirectory}/Pipeline_Log; done

ls *Bowtie3nity*.e* | while read file; do echo "" >> ${MainDirectory}/StepSummaries/SummaryOfReadSupport; echo "$file" | cut -d'.' -f1 >> ${MainDirectory}/StepSummaries/SummaryOfReadSupport ; grep -B14 "overall alignment rate" "$file" >> ${MainDirectory}/StepSummaries/SummaryOfReadSupport; echo "" | cut -d'.' -f1 >> ${MainDirectory}/StepSummaries/SummaryOfReadSupport; done

############################################
##        RETRIEVE FRAGMENT SIZES         ##
############################################

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do grep MEDIAN ${MainDirectory}/Bowtie_${id}_${project}/InsertSize_${id}_${project}.txt | tr '\t' '\n' >> ${MainDirectory}/Bowtie_${id}_${project}/colA; grep -A2 MEDIAN ${MainDirectory}/Bowtie_${id}_${project}/InsertSize_${id}_${project}.txt | grep -v MEDIAN | head -n1 | tr '\t' '\n' >> ${MainDirectory}/Bowtie_${id}_${project}/colB; grep -A2 MEDIAN ${MainDirectory}/Bowtie_${id}_${project}/InsertSize_${id}_${project}.txt | grep -v MEDIAN | tail -n1 | tr '\t' '\n' >> ${MainDirectory}/Bowtie_${id}_${project}/colC; paste ${MainDirectory}/Bowtie_${id}_${project}/colA ${MainDirectory}/Bowtie_${id}_${project}/colB ${MainDirectory}/Bowtie_${id}_${project}/colC | head -n9 >> ${MainDirectory}/Bowtie_${id}_${project}/InsertSize_Summary_${id}_${project}; rm -rf ${MainDirectory}/Bowtie_${id}_${project}/col*; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "##########  Picard Metrics For ${id}_${project}  ##########" >> ${MainDirectory}/StepSummaries/SummaryOfInsertSize; cat ${MainDirectory}/Bowtie_${id}_${project}/InsertSize_Summary_${id}_${project} >> ${MainDirectory}/StepSummaries/SummaryOfInsertSize; echo "" >> ${MainDirectory}/StepSummaries/SummaryOfInsertSize ; done; done


############################################
##          MOVE FILES TO OUTPUT          ##
############################################
mv Bowtie3nity_*.pbs ${MainDirectory}/PBSin
mv *Bowtie3nity_*.e* ${MainDirectory}/PBSout
mv *Bowtie3nity_*.o* ${MainDirectory}/PBSout
cd $MainDirectory

############################################
##             STEP  FINISHED             ##
############################################


