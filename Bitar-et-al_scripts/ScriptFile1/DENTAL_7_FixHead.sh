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

step="FIXING READ HEADERS"

echo "This script will run STEP 7 of the DENTAL pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Fixing the headers of Trimmed reads"

#  PIPELINE:
#  --------------------------------------------------------------
#   1 FASTQC PRE TRIM (quality-control before trimming)
#   2 RCORRECTOR (read correction)
#   3 FURpe (filtering uncorrectable reads)
#   4 TRIMMOMATIC (trimming reads)
#   5 FASTQC POST TRIM (quality-control after trimming)
#   6 FILTERING RIBOSOMAL READS (filter rRNA with BBduk)
#-->7 FIXING READ HEADERS (in-house script to fix formating)
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

#   7 FIXING READ HEADERS (in-house script to fix formating)
############################################
##             REASONS FOR FIX            ##
############################################

# Reason for this step: Trimmomatic seems to change the naming scheme of reads.
# Original is "/1" for left read and "/2" for right reads respectively.
# Trinity needs the naming to be kept as expected.
# Source: https://www.biostars.org/p/141602/

# WARNING!
# Single End reads that survive trimming without a pair are discarded in the current pipeline version.

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

# Create output directory:
mkdir $MainDirectory/ReadyReads

############################################
##            WRITE TO LOG FILE           ##
############################################

echo "Fixed headers directory" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/ReadyReads" >> ${MainDirectory}/Pipeline_Log
echo "What: Results for formating of read headers performed by in-house script" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) if needed:

# Memory required (in GB):
memreq=6

# CPUs required (integer):
cpureq=2

############################################
##      CREATE PBS FILES FOR FIXHEAD      ##
############################################

# PBS files for this step are named FixHead_FILE.pbs
# Cluster job names for this step are named FixHead_FILE

# List all fastq files and create correspondent PBS files:
cut -f1  ListOfFiles | while read file; do cp ${MainDirectory}/Header.pbs FixHead_${file}.pbs; done

# Populate PBS files (name of job):
cut -f1  ListOfFiles | while read file; do echo "#PBS -N FixHead_${file}" >> FixHead_${file}.pbs; done

# Populate PBS files (number of threads):
cut -f1  ListOfFiles | while read file; do echo "#PBS -r n" >> FixHead_${file}.pbs; done

# Populate PBS files (time and resources for running the job):
cut -f1  ListOfFiles | while read file; do echo "#PBS -l mem=${memreq}GB,walltime=24:00:00,ncpus=${cpureq}" >> FixHead_${file}.pbs; done

# Populate PBS files (send message when job is finished or aborted):
cut -f1  ListOfFiles | while read file; do echo "#PBS -m ae" >> FixHead_${file}.pbs; done

# Populate PBS files (e-mail address for correspondence):
#cut -f1  ListOfFiles | while read file; do echo "#PBS -M $email" >> FixHead_${file}.pbs; done

# Populate PBS files (empty line):
cut -f1  ListOfFiles | while read file; do echo "" >> FixHead_${file}.pbs; done

############################################
##        COMMAND LINE FOR FIXHEAD        ##
############################################
# Populate PBS files (actual command lines):

# Paired and Trimmed reads:
cut -f1 ListOfFiles | while read file; do f1="${MainDirectory}/RiboDepleted/Trimmed_${file}_R1.fq"; echo "cat ${f1} | awk '{ if (NR%4==1) { print"' $1substr($2,2)"/1" }'" else { print } }' | sed 's/_1://g' >> ${MainDirectory}/ReadyReads/CleanUNcorrected_${file}_R1.fq; gzip ${MainDirectory}/ReadyReads/CleanUNcorrected_${file}_R1.fq; gzip ${f1}" >> FixHead_${file}.pbs ; done

cut -f1 ListOfFiles | while read file; do f2="${MainDirectory}/RiboDepleted/Trimmed_${file}_R2.fq"; echo "cat ${f2} | awk '{ if (NR%4==1) { print"' $1substr($2,2)"/2" }'" else { print } }' | sed 's/_2://g' >> ${MainDirectory}/ReadyReads/CleanUNcorrected_${file}_R2.fq; gzip ${MainDirectory}/ReadyReads/CleanUNcorrected_${file}_R2.fq; gzip ${f2}" >> FixHead_${file}.pbs ; done

# Paired, Corrected and Trimmed reads:
cut -f1 ListOfFiles | while read file; do f1="${MainDirectory}/RiboDepleted/CorTrimmed_${file}_R1.cor.fq"; echo "cat ${f1} | awk '{ if (NR%4==1) { print"' $1substr($2,2)"/1" }'" else { print } }' | sed 's/_1://g' >> ${MainDirectory}/ReadyReads/CleanCorrected_${file}_R1.cor.fq; gzip ${MainDirectory}/ReadyReads/CleanCorrected_${file}_R1.cor.fq; gzip ${f1}" >> FixHead_${file}.pbs ; done

cut -f1 ListOfFiles | while read file; do f2="${MainDirectory}/RiboDepleted/CorTrimmed_${file}_R2.cor.fq"; echo "cat ${f2} | awk '{ if (NR%4==1) { print"' $1substr($2,2)"/2" }'" else { print } }' | sed 's/_2://g' >> ${MainDirectory}/ReadyReads/CleanCorrected_${file}_R2.cor.fq; gzip ${MainDirectory}/ReadyReads/CleanCorrected_${file}_R2.cor.fq; gzip ${f2}" >> FixHead_${file}.pbs ; done


############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "===> Read header formating" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- command lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n4 FixHead*.pbs >> $MainDirectory/Pipeline_Log

############################################
##          SUBMIT FIXHEAD JOBS           ##
############################################
# Populate Log file (header):
echo "Job Submission IDs:" >> ${MainDirectory}/PBSsub/SubJobs_FixHead

# Submit jobs and populate Log file (submission ID numbers):
ls FixHead_*.pbs | while read line; do echo "$line" >> ${MainDirectory}/PBSsub/SubJobs_FixHead ; qsub "$line"; done >> ${MainDirectory}/PBSsub/SubJobs_FixHead

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user" | grep FixHead
qstat | grep -c FixHead

############################################
##                   END                  ##
############################################

