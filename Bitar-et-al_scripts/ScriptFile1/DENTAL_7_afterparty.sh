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
##          MOVE FILES TO OUTPUT          ##
############################################
cd $MainDirectory
mv FixHead*.pbs ${MainDirectory}/PBSin
mv *FixHead*.e* ${MainDirectory}/PBSout
mv *FixHead*.o* ${MainDirectory}/PBSout

############################################
##             STEP  FINISHED             ##
############################################

############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""

