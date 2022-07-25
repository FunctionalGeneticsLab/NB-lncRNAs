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
ls *RSeQC*.o* | while read line; do echo "$line" >> ${MainDirectory}/StepSummaries/SummaryOfStrandness ; cat "$line" >> ${MainDirectory}/StepSummaries/SummaryOfStrandness; echo "******************************" >> ${MainDirectory}/StepSummaries/SummaryOfStrandness; done

############################################
##          MOVE FILES TO OUTPUT          ##
############################################
mv RSeQC*.pbs ${MainDirectory}/PBSin
mv *RSeQC*.e* ${MainDirectory}/PBSout
mv *RSeQC*.o* ${MainDirectory}/PBSout
cd $MainDirectory

############################################
##             STEP  FINISHED             ##
############################################

############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""

