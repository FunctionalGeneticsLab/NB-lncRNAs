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



############################################
##          MOVE FILES TO OUTPUT          ##
############################################
cd $MainDirectory
mv *Trinity*.pbs $MainDirectory/PBSin
mv *Trinity*.e* $MainDirectory/PBSout
mv *Trinity*.o* $MainDirectory/PBSout

# Centralise output assemblies to the main Trinity folder:
mkdir ${MainDirectory}/Trinity
mv ${MainDirectory}/*.Trinity.fasta ${MainDirectory}/Trinity
mv ${MainDirectory}/*.Trinity.fasta.gene_trans_map ${MainDirectory}/Trinity


############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""

