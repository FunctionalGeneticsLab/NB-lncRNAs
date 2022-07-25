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
##          MERGE OUTPUT FOLDERS          ##
############################################

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do mkdir ${MainDirectory}/Busco_${id}_${project}; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do mv ${MainDirectory}/run_Busco*_${id}_${project} ${MainDirectory}/Busco_${id}_${project}; mv ${MainDirectory}/BuscoTemp*_${id}_${project} ${MainDirectory}/Busco_${id}_${project}; done; done

############################################
##              MAKE SUMMARY              ##
############################################
# See BUSCO results:
cd $MainDirectory

ls run_BuscoMammalia_* | grep ':' | tr -d ':' | cut -d'_' -f2- | while read folder; do grep "Complete" run_${folder}/short_summary_${folder}.txt | head -n1 ; done
ls run_BuscoEukaryota_* | grep ':' | tr -d ':' | cut -d'_' -f2- | while read folder; do grep "Complete" run_${folder}/short_summary_${folder}.txt | head -n1 ; done

# Summarize BUSCO results:
echo "** Busco Results **" >> ${MainDirectory}/Pipeline_Log

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> ${MainDirectory}/StepSummaries/SummaryOfBusco; echo "##### Eukaryota BUSCO: ${id}_${project} #####" >> ${MainDirectory}/StepSummaries/SummaryOfBusco ; tail ${MainDirectory}/Busco_${id}_${project}/run_BuscoEukaryota_${id}_${project}/short_summary_BuscoEukaryota_${id}_${project}.txt >> ${MainDirectory}/StepSummaries/SummaryOfBusco; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> ${MainDirectory}/StepSummaries/SummaryOfBusco; echo "##### Mammalia BUSCO: ${id}_${project} #####" >> ${MainDirectory}/StepSummaries/SummaryOfBusco ; tail ${MainDirectory}/Busco_${id}_${project}/run_BuscoMammalia_${id}_${project}/short_summary_BuscoMammalia_${id}_${project}.txt >> ${MainDirectory}/StepSummaries/SummaryOfBusco; done; done

cat ${MainDirectory}/StepSummaries/SummaryOfBusco >> ${MainDirectory}/Pipeline_Log

############################################
##          MOVE FILES TO OUTPUT          ##
############################################
mv Busco_*.pbs ${MainDirectory}/PBSin
mv *Busco_*.e* ${MainDirectory}/PBSout
mv *Busco_*.o* ${MainDirectory}/PBSout

rm -rf BuscoTemp*

cd $MainDirectory
############################################
##             STEP  FINISHED             ##
############################################


