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
##            DECLARE  MODULES            ##
############################################

############################################
##              MAKE SUMMARY              ##
############################################
# Summarize Genome Alignment results:
cd $MainDirectory
ls *Bowtie*.e* | while read line; do echo "$line" >> ${MainDirectory}/StepSummaries/SummaryOfGenomeAlignment ; echo "" >> ${MainDirectory}/StepSummaries/SummaryOfGenomeAlignment ; head -n5 "$line" >> ${MainDirectory}/StepSummaries/SummaryOfGenomeAlignment; echo "******************************" >> ${MainDirectory}/StepSummaries/SummaryOfGenomeAlignment; grep  "overall alignment rate" "$line" >> ${MainDirectory}/StepSummaries/SummaryOfGenomeAlignment;  echo "******************************" >> ${MainDirectory}/StepSummaries/SummaryOfGenomeAlignment; tail "$line" | grep -v "overall alignment rate" >> ${MainDirectory}/StepSummaries/SummaryOfGenomeAlignment; echo "" >> ${MainDirectory}/StepSummaries/SummaryOfGenomeAlignment; done

############################################
##          MOVE FILES TO OUTPUT          ##
############################################
mv Bowtie*.pbs ${MainDirectory}/PBSin
mv *Bowtie*.e* ${MainDirectory}/PBSout
mv *Bowtie*.o* ${MainDirectory}/PBSout
cd $MainDirectory

############################################
##             STEP  FINISHED             ##
############################################

############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""

