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
# Summarize Trimmomatic results:
cd ${MainDirectory}

ls *Trimmomatic*.e* | while read line; do echo "$line" >> ${MainDirectory}/StepSummaries/SummaryOfBBduk ; grep -A10 "^Input is being processed as paired" "$line" | grep ":" >> ${MainDirectory}/StepSummaries/SummaryOfBBduk; echo "" >> ${MainDirectory}/StepSummaries/SummaryOfBBduk; done

ls *Trimmomatic*.e* | while read line; do echo "$line" >> ${MainDirectory}/StepSummaries/SummaryOfTrimming ; grep -A2 "^ILLUMINACLIP:" "$line" >> ${MainDirectory}/StepSummaries/SummaryOfTrimming; echo "" >> ${MainDirectory}/StepSummaries/SummaryOfTrimming; done

# Additionally, one may use this command to check that all jobs finished successfully:
# grep -c "TrimmomaticPE: Completed successfully" ${MainDirectory}/StepSummaries/SummaryOfTrimming

# Calculate statistics for filtering steps (the Statistics.sh script must be available):
# In-house Statistics script:
echo "********************" >> ${MainDirectory}/StepSummaries/SummaryOfTrimming
echo "Percentage of Surviving Pairs:" >> ${MainDirectory}/StepSummaries/SummaryOfTrimming
StatScript="/working/lab_julietF/mainaB/Statistics.sh"
grep "Input Read Pairs: " ${MainDirectory}/StepSummaries/SummaryOfTrimming | cut -d')' -f1 | cut -d'(' -f2 | tr -d '%' | ${StatScript} >> ${MainDirectory}/StepSummaries/SummaryOfTrimming


############################################
##          MOVE FILES TO OUTPUT          ##
############################################
cd $MainDirectory
mv Trimmomatic*.pbs ${MainDirectory}/PBSin
mv *Trimmomatic*.e* ${MainDirectory}/PBSout
mv *Trimmomatic*.o* ${MainDirectory}/PBSout

############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""

