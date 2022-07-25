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

step="FASTQC POST TRIMMING"

echo "This script will run STEP 5 of the DENTAL pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Perform QC of reads after trimming"

#  PIPELINE:
#  --------------------------------------------------------------
#   1 FASTQC PRE TRIM (quality-control before trimming)
#   2 RCORRECTOR (read correction)
#   3 FURpe (filtering uncorrectable reads)
#   4 TRIMMOMATIC (trimming reads)
#-->5 FASTQC POST TRIM (quality-control after trimming)
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

#   5 FASTQC POST TRIM (quality-control after trimming)
############################################
##             SOURCE OF TOOLS            ##
############################################

# FastQC was obtained from:
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc

############################################
##            REQUIRED MODULES            ##
############################################
python=python/2.7.10
python3=python/3.6.1
fastqc=fastqc/0.11.8

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

# Create FastQC-PostTrim directory:
mkdir $MainDirectory/FastQCPostTrim

############################################
##            WRITE TO LOG FILE           ##
############################################

echo "FastQC-PostTrim directory" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/FastQCPostTrim" >> ${MainDirectory}/Pipeline_Log
echo "What: Results for quality control performed by FastQC on trimmed reads" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) if needed:

# Memory required (in GB):
memreq=10

# CPUs required (integer):
cpureq=4

############################################
##      CREATE PBS FILES FOR FASTQC       ##
############################################

# List all fastq files and create correspondent PBS files:
cut -f1 ListOfFiles | while read file; do cp $MainDirectory/Header.pbs postTrimFastQC_${file}.pbs; done

# Populate PBS files (name of job):
cut -f1 ListOfFiles | while read file; do echo "#PBS -N ${file}_1FastQC" >> postTrimFastQC_${file}.pbs; done

# Populate PBS files (number of threads):
cut -f1 ListOfFiles | while read file; do echo "#PBS -r n" >> postTrimFastQC_${file}.pbs; done

# Populate PBS files (time and resources for running the job):
cut -f1 ListOfFiles | while read file; do echo "#PBS -l mem=${memreq}GB,walltime=08:00:00,ncpus=${cpureq}" >> postTrimFastQC_${file}.pbs; done

# Populate PBS files (send message when job is finished or aborted):
cut -f1 ListOfFiles | while read file; do echo "#PBS -m ae" >> postTrimFastQC_${file}.pbs; done

# Populate PBS files (e-mail address for correspondence):
cut -f1 ListOfFiles | while read file; do echo "#PBS -M $email" >> postTrimFastQC_${file}.pbs; done

# Populate PBS files (empty line):
cut -f1 ListOfFiles | while read file; do echo "" >> postTrimFastQC_${file}.pbs; done

# Populate PBS files (load module(s)):
cut -f1 ListOfFiles | while read file; do echo "module load $fastqc" >> postTrimFastQC_${file}.pbs; done

# Populate PBS files (empty line):
cut -f1 ListOfFiles | while read file; do echo "" >> postTrimFastQC_${file}.pbs; done

############################################
##         COMMAND LINE FOR FASTQC        ##
############################################
# Populate PBS files (actual command line):
cut -f1 ListOfFiles | while read file; do echo "fastqc -t $cpureq --outdir ${MainDirectory}/FastQCPostTrim $MainDirectory/Trimmomatic/CorTrimmed_${file}_R1.fq.gz $MainDirectory/Trimmomatic/CorTrimmed_${file}_R2.fq.gz" >> postTrimFastQC_${file}.pbs ; done

############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "===> FastQC Post Trimming" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- command lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n1 postTrimFastQC_*.pbs >> $MainDirectory/Pipeline_Log

############################################
##            SUBMIT FASTQC JOBS          ##
############################################
# Populate Log file (header):
echo "-- FastQC (PostTrim) Job Submission IDs --" >> $MainDirectory/PBSsub/SubJobs_postTrimFastQC

# Submit jobs and populate Log file (submission ID numbers):
ls postTrimFastQC_*.pbs | while read line; do echo "$line" >> $MainDirectory/PBSsub/SubJobs_postTrimFastQC ; qsub "$line"; done >> $MainDirectory/PBSsub/SubJobs_postTrimFastQC

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user" | grep Fast
qstat | grep -c Fast

############################################
##                   END                  ##
############################################

