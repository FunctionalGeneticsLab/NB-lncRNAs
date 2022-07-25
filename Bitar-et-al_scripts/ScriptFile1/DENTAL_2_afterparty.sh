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
# Summarize Rcorrector results:
cd ${MainDirectory}
grep -B1 -Pe "^\tCorrected " Rcorrector*.e* | tr '-' '\t' | tr -d ':' >> ${MainDirectory}/StepSummaries/SummaryOfRcorrector

############################################
##          MOVE FILES TO OUTPUT          ##
############################################
mv Rcorrector*.pbs ${MainDirectory}/PBSin
mv *Rcorrector*.e* ${MainDirectory}/PBSout
mv *Rcorrector*.o* ${MainDirectory}/PBSout
rm -rf ${MainDirectory}/TMP_Rcorrector
cd $MainDirectory

############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""

