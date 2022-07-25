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
cut -f1 ListOfFiles | while read file; do echo "$file" >> ${MainDirectory}/StepSummaries/SummaryOfFURpe ; cat ${MainDirectory}/FURpe/rmunfixable_${file}.log >> ${MainDirectory}/StepSummaries/SummaryOfFURpe; echo "**********" >> ${MainDirectory}/StepSummaries/SummaryOfFURpe; echo "" >> ${MainDirectory}/FURpe/SummaryUnfixable; done

cd ${MainDirectory}/FURpe
grep "^S" SummaryUnfixable >> Col1; grep "^total" SummaryUnfixable | cut -d':' -f2 >> Col2; grep "^removed" SummaryUnfixable | cut -d':' -f2 >> Col3; grep "^retained" SummaryUnfixable | cut -d':' -f2 >> Col4; grep "^pairs" SummaryUnfixable | cut -d':' -f2 >> Col5; grep "^both" SummaryUnfixable | cut -d':' -f2 >> Col6; echo "SampleName" >> ColA; grep "^total" SummaryUnfixable | cut -d':' -f1 | uniq >> ColB; grep "^removed" SummaryUnfixable | cut -d':' -f1 | uniq >> ColC; grep "^retained" SummaryUnfixable | cut -d':' -f1 | uniq >> ColD; grep "^pairs" SummaryUnfixable | cut -d':' -f1 | uniq >> ColE; grep "^both" SummaryUnfixable | cut -d':' -f1 | uniq >> ColF; paste ColA ColB ColC ColD ColE ColF | tr ' ' '_' >> ${MainDirectory}/StepSummaries/SummaryOfFURpe ; paste Col1 Col2 Col3 Col4 Col5 Col6 >> ${MainDirectory}/StepSummaries/SummaryOfFURpe; rm -rf Col*


############################################
##          MOVE FILES TO OUTPUT          ##
############################################
cd $MainDirectory
mv FURpe*.pbs ${MainDirectory}/PBSin
mv *FURpe*.e* ${MainDirectory}/PBSout
mv *FURpe*.o* ${MainDirectory}/PBSout

############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""

