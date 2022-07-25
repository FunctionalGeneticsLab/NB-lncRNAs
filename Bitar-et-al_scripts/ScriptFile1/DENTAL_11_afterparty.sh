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

# In-house Statistics script:
StatScript="/working/lab_julietF/mainaB/Scripts/Statistics.sh"

# Spearman correlation:
cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do max=`grep -v "#" ${MainDirectory}/DeepTools/SpearmanCorrelation_${project}.matrix | head -n1 | tr '\t' '\n' | wc -l`; seq 2 ${max} | while read n; do grep -v "#" ${MainDirectory}/DeepTools/SpearmanCorrelation_${project}.matrix | cut -f${n} | head -n1 | while read sample; do grep "$sample" ${MainDirectory}/DeepTools/SpearmanCorrelation_${project}.matrix | cut -f1,${n} >> ${MainDirectory}/DeepTools/SpearmanCorrelation_${project}.samples; done; done; done

cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${project}"; cut -f2 ${MainDirectory}/DeepTools/SpearmanCorrelation_${project}.samples | grep -v "_\|1.0000" | ${StatScript} ; done

cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${project}"; cut -f2- DeepTools/SpearmanCorrelation_${project}.matrix | grep -v "^#\|_" | tr '\t' '\n' | grep -v "^1.0000" | ${StatScript}; done


# Pearson correlation:
cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do max=`grep -v "#" ${MainDirectory}/DeepTools/PearsonCorrelation_${project}.matrix | head -n1 | tr '\t' '\n' | wc -l`; seq 2 ${max} | while read n; do grep -v "#" ${MainDirectory}/DeepTools/PearsonCorrelation_${project}.matrix | cut -f${n} | head -n1 | while read sample; do grep "$sample" ${MainDirectory}/DeepTools/PearsonCorrelation_${project}.matrix | cut -f1,${n} >> ${MainDirectory}/DeepTools/PearsonCorrelation_${project}.samples; done; done; done

cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${project}"; cut -f2 ${MainDirectory}/DeepTools/PearsonCorrelation_${project}.samples | grep -v "_\|1.0000" | ${StatScript} ; done

cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${project}"; cut -f2- DeepTools/PearsonCorrelation_${project}.matrix | grep -v "^#\|_" | tr '\t' '\n' | grep -v "^1.0000" | ${StatScript} ; done


############################################
##          MOVE FILES TO OUTPUT          ##
############################################
cd $MainDirectory
mv deeptools*.pbs $MainDirectory/PBSin
mv deeptools*.e* $MainDirectory/PBSout
mv deeptools*.o* $MainDirectory/PBSout


############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""

