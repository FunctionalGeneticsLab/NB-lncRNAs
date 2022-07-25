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
# Summarize Ribodepletion results:
ls *Ribosomal*.e* | while read line; do echo "$line" >> ${MainDirectory}/StepSummaries/SummaryOfRiboDepletion ; echo "(1) Paired Trimmed - (2) Paired Corrected and Trimmed - (3) Unpaired Trimmed - (4) Unpaired Corrected and Trimmed" >> ${MainDirectory}/StepSummaries/SummaryOfRiboDepletion ; echo "" >> ${MainDirectory}/StepSummaries/SummaryOfRiboDepletion ; grep -A1 "^Result: \|^Contaminants: " "$line" | grep -v "^Total Removed:\|^--" >> ${MainDirectory}/StepSummaries/SummaryOfRiboDepletion; echo "******************************" >> ${MainDirectory}/StepSummaries/SummaryOfRiboDepletion; done

############################################
##          MOVE FILES TO OUTPUT          ##
############################################
cd $MainDirectory
mv Ribosomal*.pbs ${MainDirectory}/PBSin
mv *Ribosomal*.e* ${MainDirectory}/PBSout
mv *Ribosomal*.o* ${MainDirectory}/PBSout

############################################
##             STEP  FINISHED             ##
############################################

############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""

