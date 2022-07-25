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
##             PIPELINE STEPS             ##
############################################

step="READ CORRECTION"

echo "This script will run STEP 2 of the DENTAL pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Correct reads using Rcorrector"

#  PIPELINE:
#  --------------------------------------------------------------
#   1 FASTQC PRE TRIM (quality-control before trimming)
#   2 RCORRECTOR (read correction)
#   3 FILTER OUT UNCORRECTABLE READS (removed with FUR-pe)
#   4 TRIMMOMATIC (trimming reads)
#   5 FASTQC POST TRIM (quality-control after trimming)
#   6 FILTERING RIBOSOMAL READS (filter rRNA with BBduk)
#   7 FIXING READ HEADERS (in-house script to fix formating)
#   8 MAP TO GENOME (using Bowtie to flag contamination)
#   9 INFER EXPERIMENTAL DESIGN (infer design using RSeQC)
#   10 FASTQC PRE TRINITY (final quality-control before Trinity)
#  --------------------------------------------------------------
#   11 ASSESS REPLICATES CORRELATION
#   12 TRINITY
#   13 GENERAL QUALITY STATISTICS
#   14 TRANSCRIPTOME READ REPRESENTATION
#   15 BUSCO (Benchmarking Universal Single-Copy Orthologs)
#   16 TransRate (Analyse de novo assembly)
#   17 Detonate (DE novo TranscriptOme rNa-seq Assembly with or without the Truth Evaluation)
#   18 ALIGN TRANSCRIPTS TO GENOME AND ASSESS SPLICING
#   19 FILTER OUT LOW SUPPORT TRANSCRIPTS
#   20 MERGE TRANSCRIPTS TO SUPER TRANSCRIPTS
#  --------------------------------------------------------------
#   21 RETRIEVE ANNOTATION
#   22 PREDICT CODING POTENTIAL
#   23 Identify ncRNAs with FEELnc
#   24 COMPARE ALL TOOLS DEFINE lncRNAs
#   25 METRICS OF FINAL lncRNAs
#  --------------------------------------------------------------
#   26 COMPARE TRINITY ASSEMBLIES
#   27 PREPARE FOR DE ANALYSES
#   28 BUILD SALMON INDICES
#   29 DE ANALYSES WITH SALMON
#   30 PERFORM DE ANALYSES WITH EDGER
#  --------------------------------------------------------------


############################################
##            EXAMPLES OF DATA            ##
############################################

# Data should be named as:
# ID1_ID2_R1.fastq.gz (forward/reverse)
# ID1_ID2_R2.fastq.gz (reverse/forward)

############################################
##              SET USER VARs             ##
############################################
# Edit the following file(s) to reflect your environment

# Email to which messages will be sent when each PBS job is initiated/concluded/terminated
email=Maina.Bitar@qimrberghofer.edu.au

# User name within your cluster environment
user=`whoami`

############################################
##         SET MAIN WORK DIRECTORY        ##
############################################
# Define directory:

if [ "$#" -gt 0 ]; then MainDirectory=$1; fi
if [ "$#" -eq 0 ]; then MainDirectory=`pwd`; fi

echo ""; echo "WARNING: The main directory for this run was set to ${MainDirectory}."; echo ""

############################################
##        CREATE OUTPUT DIRECTORIES       ##
############################################

# Go inside Main Directory (containing ONLY original fastq files)
cd $MainDirectory

# Create Main Output directories:
mkdir ${MainDirectory}/PBSin
mkdir ${MainDirectory}/PBSout
mkdir ${MainDirectory}/PBSsub
mkdir ${MainDirectory}/StepSummaries


############################################
##          CREATE FILELIST FILE          ##
############################################
# Go inside Main Directory (containing ONLY original fastq files):
cd $MainDirectory

# Assuming files are named following the example, run:

ls *fq.gz | cut -d'_' -f1-2 | sort | uniq >> ListOfFiles

#ls *fastq.gz | cut -d'_' -f1-2 | sort | uniq >> ListOfFiles

############################################
##           CREATE HEADER FILE           ##
############################################
# Populate header file:
echo "##########################################################################" >> ${MainDirectory}/Header.pbs
echo "#" >> ${MainDirectory}/Header.pbs
echo "#  Script:    Script for DENTAL - DE Novo Transcriptome Assembly Line" >> ${MainDirectory}/Header.pbs
echo "#  Author:    Maina Bitar" >> ${MainDirectory}/Header.pbs
echo "#  Created:   2022 at QIMR Berghofer [Brisbane, Australia]" >> ${MainDirectory}/Header.pbs
echo "#  Email:     Maina.Bitar@qimrberghofer.edu.au" >> ${MainDirectory}/Header.pbs
echo "#" >> ${MainDirectory}/Header.pbs
echo "##########################################################################" >> ${MainDirectory}/Header.pbs

############################################
##             CREATE LOG FILE            ##
############################################
# Populate Log file:

echo "***-- WELCOME --***" >> ${MainDirectory}/Pipeline_Log
echo "To DENTAL - DE Novo Transcriptome Assembly Line" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

echo -n "It is now exactly " >> ${MainDirectory}/Pipeline_Log
date >> ${MainDirectory}/Pipeline_Log
echo "And we are just getting started..." >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

echo "*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***" >> ${MainDirectory}/Pipeline_Log
echo "We created folders for your PBS inputs and outputs:" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

echo "PBS input directory" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/PBSin" >> ${MainDirectory}/Pipeline_Log
echo "What: Stores all PBS job submission scripts, therefore all actual command lines" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

echo "PBS output directory" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/PBSout" >> ${MainDirectory}/Pipeline_Log
echo "What: Stores all PBS job submission errors and output messages" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

echo "*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***" >> ${MainDirectory}/Pipeline_Log
echo "We will create new folders for the outputs of PART 1:" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

echo "FastQC-PreTrim directory" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/FastQCPreTrim" >> ${MainDirectory}/Pipeline_Log
echo "What: Results for quality control performed by FastQC on pre-trimming reads" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

echo "Rcorrector directory" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/Rcorrector" >> ${MainDirectory}/Pipeline_Log
echo "What: Results of correcting reads with Rcorrector" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

echo "FURpe directory" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/FURpe" >> ${MainDirectory}/Pipeline_Log
echo "What: Results of filtering uncorrectable reads using FURpe" >> ${MainDirectory}/Pipeline_Log
echo "FilterUncorrectabledPEfastq.py" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

echo "FastQCPostTrim directory" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/FastQCPostTrim" >> ${MainDirectory}/Pipeline_Log
echo "What: Results for quality control performed by FastQC on post-trimming reads" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

echo "Trimmomatic directory" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/Trimmomatic" >> ${MainDirectory}/Pipeline_Log
echo "What: Results for trimming performed by Trimmomatic on all reads" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

echo "FastQCPostTrim directory" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/FastQCPostTrim" >> ${MainDirectory}/Pipeline_Log
echo "What: Results for quality control performed by FastQC on post-trimming reads" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

echo "*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***" >> ${MainDirectory}/Pipeline_Log


############################################
##                   END                  ##
############################################

echo ""; echo "---> FINISHED"; echo ""
