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
python=python/2.7.10
python3=python/3.5.5
fastqc=fastqc/0.11.8

############################################
##          MOVE FILES TO OUTPUT          ##
############################################
cd $MainDirectory
mv postTrimFastQC*.pbs $MainDirectory/PBSin
mv *FastQC*.e* $MainDirectory/PBSout
mv *FastQC*.o* $MainDirectory/PBSout

############################################
##               RUN MULTIQC              ##
############################################
cd $MainDirectory/FastQCPostTrim

echo "Running MultiQC to consolidate results for these files:"
ls $MainDirectory/FastQCPostTrim/*html | cut -d'.' -f1

module load $python3
multiqc ./

############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""

