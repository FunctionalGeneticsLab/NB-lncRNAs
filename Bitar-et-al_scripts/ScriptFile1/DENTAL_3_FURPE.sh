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

step="FILTER UNCORRECTED READS"

echo "This script will run STEP 3 of the DENTAL pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Filter out uncorrectable reads with FUR-pe"

#  PIPELINE:
#  --------------------------------------------------------------
#   1 FASTQC PRE TRIM (quality-control before trimming)
#   2 RCORRECTOR (read correction)
#-—>3 FURpe (filtering uncorrectable reads)
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

#   3 FURpe (filtering uncorrectable reads)
############################################
##             SOURCE OF TOOLS            ##
############################################

# FURpe (FilterUncorrectabledPEfastq) was written by Adam Freedman and was obtained from:
# https://github.com/harvardinformatics/TranscriptomeAssemblyTools/blob/master/FilterUncorrectabledPEfastq.py

# The original script was converted to python3 format using 2to3 and both versions can be used.

############################################
##             TOOLS DIRECTORY            ##
############################################
# User-installed tools
# Local directory where all tools below are installed:
toolsdir="/working/lab_julietF/mainaB/Tools"

############################################
##            REQUIRED MODULES            ##
############################################
furpe="${toolsdir}/FilterUncorrectabledPEfastq_python2.py"

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

# Create FURpe directory:
mkdir ${MainDirectory}/FURpe

############################################
##            WRITE TO LOG FILE           ##
############################################

echo "FURpe directory" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/FURpe" >> ${MainDirectory}/Pipeline_Log
echo "What: Results of filtering uncorrectable reads using FURpe" >> ${MainDirectory}/Pipeline_Log
echo "FilterUncorrectabledPEfastq.py" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) to reflect your requirements:

# Memory required (in GB):
memreq=2

# CPUs required (integer):
cpureq=2


############################################
##      CREATE PBS FILES FOR FUR-pe       ##
############################################
# PBS files for this step are named FURpe_FILE.pbs
# Cluster job names for this step are named FURpe_FILE

# List all fastq files and create correspondent PBS files:
cut -f1 ListOfFiles | while read file; do cp ${MainDirectory}/Header.pbs FURpe_${file}.pbs; done

# Populate PBS files (name of job):
cut -f1 ListOfFiles | while read file; do echo "#PBS -N FURpe_${file}" >> FURpe_${file}.pbs; done

# Populate PBS files (number of threads):
cut -f1 ListOfFiles | while read file; do echo "#PBS -r n" >> FURpe_${file}.pbs; done

# Populate PBS files (time and resources for running the job):
cut -f1 ListOfFiles | while read file; do echo "#PBS -l mem=${memreq}GB,walltime=08:00:00,ncpus=${cpureq}" >> FURpe_${file}.pbs; done

# Populate PBS files (send message when job is finished or aborted):
cut -f1 ListOfFiles | while read file; do echo "#PBS -m ae" >> FURpe_${file}.pbs; done

# Populate PBS files (e-mail address for correspondence):
cut -f1 ListOfFiles | while read file; do echo "#PBS -M $email" >> FURpe_${file}.pbs; done

# Populate PBS files (change directory):
cut -f1 ListOfFiles | while read file; do echo "cd ${MainDirectory}/FURpe" >> FURpe_${file}.pbs; done

# Populate PBS files (empty line):
cut -f1 ListOfFiles | while read file; do echo "" >> FURpe_${file}.pbs; done

# Populate PBS files (load module(s)):
cut -f1 ListOfFiles | while read file; do echo "module load $python2" >> FURpe_${file}.pbs; done

# Populate PBS files (empty line):
cut -f1 ListOfFiles | while read file; do echo "" >> FURpe_${file}.pbs; done

############################################
##     COMMAND LINE FOR FUR-pe filter     ##
############################################
# Populate PBS files (actual command line):
cut -f1 ListOfFiles | while read file; do f1="${MainDirectory}/Rcorrector/${file}/${file}_R1.cor.fq.gz"; f2="${MainDirectory}/Rcorrector/${file}/${file}_R2.cor.fq.gz"; echo "python $furpe -1 ${f1} -2 ${f2} -s ${file}" >> FURpe_${file}.pbs ; done

############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (empty line):
echo "" >> ${MainDirectory}/Pipeline_Log

# Populate Log file (header):
echo "** FURpe Command Lines **" >> ${MainDirectory}/Pipeline_Log

# Populate Log file (command lines used):
tail -n1 FURpe_*.pbs >> ${MainDirectory}/Pipeline_Log

# Populate Log file (empty line):
echo "" >> ${MainDirectory}/Pipeline_Log

############################################
##           SUBMIT FUR-pe JOBS           ##
############################################

# Populate Log file (header):
echo "FURpe Job Submission IDs:" >> ${MainDirectory}/PBSsub/SubJobs_FURpe

# Submit jobs and populate Log file (submission ID numbers):
ls FURpe_*.pbs | while read line; do echo "$line" >> ${MainDirectory}/PBSsub/SubJobs_FURpe ; qsub "$line"; done >> ${MainDirectory}/PBSsub/SubJobs_FURpe

# Inform status of jobs after submission:
qstat | grep "$user"
qstat | grep -c FURpe

############################################
##                   END                  ##
############################################
