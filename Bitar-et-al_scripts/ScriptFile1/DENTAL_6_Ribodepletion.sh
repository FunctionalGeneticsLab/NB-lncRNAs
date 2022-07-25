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

step="FILTERING RIBOSOMAL READS"

echo "This script will run STEP 6 of the DENTAL pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Filter rRNA-derived reads with BBduk"

#  PIPELINE:
#  --------------------------------------------------------------
#   1 FASTQC PRE TRIM (quality-control before trimming)
#   2 RCORRECTOR (read correction)
#   3 FURpe (filtering uncorrectable reads)
#   4 TRIMMOMATIC (trimming reads)
#   5 FASTQC POST TRIM (quality-control after trimming)
#-->6 FILTERING RIBOSOMAL READS (filter rRNA with BBduk)
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

#   6 FILTERING RIBOSOMAL READS (filter rRNA with BBduk)
############################################
##             SOURCE OF TOOLS            ##
############################################

# BBDuk is part of BBTools:
# https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/

############################################
##             TOOLS DIRECTORY            ##
############################################
# User-installed tools
# Local directory where all tools below are installed:
toolsdir="/working/lab_julietF/mainaB/Tools"

############################################
##            REQUIRED MODULES            ##
############################################
bbduk="${toolsdir}/bbmap/bbduk.sh"

############################################
##           SET REFERENCE VARs           ##
############################################
# Main directory serving as repository for references:
referencedir=/working/lab_julietF/mainaB/ReferenceGenomes

# rRNA reference:
ref=${referencedir}/Combined_human_ribosomal.fa

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

# Create BBDuk ribodepletion directory:
mkdir $MainDirectory/RiboDepleted

############################################
##            WRITE TO LOG FILE           ##
############################################

echo "BBDuk ribodepletion directory" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/RiboDepleted" >> ${MainDirectory}/Pipeline_Log
echo "What: Results for depletion of rRNA reads performed by BBDuk" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) if needed:

# Memory required (in GB):
memreq=16
memless=12

# CPUs required (integer):
cpureq=8
cpuless=5

############################################
##       CREATE PBS FILES FOR BBDuk       ##
############################################

# List all fastq files and create correspondent PBS files:
cut -f1  ListOfFiles | while read file; do cp ${MainDirectory}/Header.pbs Ribosomal_${file}.pbs; done

# Populate PBS files (name of job):
cut -f1  ListOfFiles | while read file; do echo "#PBS -N Ribosomal_${file}" >> Ribosomal_${file}.pbs; done

# Populate PBS files (number of threads):
cut -f1  ListOfFiles | while read file; do echo "#PBS -r n" >> Ribosomal_${file}.pbs; done

# Populate PBS files (time and resources for running the job):
cut -f1  ListOfFiles | while read file; do echo "#PBS -l mem=${memreq}GB,walltime=08:00:00,ncpus=${cpureq}" >> Ribosomal_${file}.pbs; done

# Populate PBS files (send message when job is finished or aborted):
cut -f1  ListOfFiles | while read file; do echo "#PBS -m ae" >> Ribosomal_${file}.pbs; done

# Populate PBS files (e-mail address for correspondence):
cut -f1  ListOfFiles | while read file; do echo "#PBS -M $email" >> Ribosomal_${file}.pbs; done

# Populate PBS files (empty line):
cut -f1  ListOfFiles | while read file; do echo "" >> Ribosomal_${file}.pbs; done

############################################
##         COMMAND LINE FOR BBDuk         ##
############################################
# Populate PBS files (actual command line):

# Paired and Trimmed reads:
cut -f1 ListOfFiles | while read file; do f1="${MainDirectory}/Trimmomatic/Trimmed_${file}_R1.fq.gz"; f2="${MainDirectory}/Trimmomatic/Trimmed_${file}_R2.fq.gz"; echo "$toolsdir/bbmap/bbduk.sh t=$cpuless -Xmx${memless}g in=${f1} in2=${f2} outm=${MainDirectory}/RiboDepleted/riboRNA_${file}_R1.fq.gz outm2=${MainDirectory}/RiboDepleted/riboRNA_${file}_R2.fq.gz out=${MainDirectory}/RiboDepleted/Trimmed_${file}_R1.fq out2=${MainDirectory}/RiboDepleted/Trimmed_${file}_R2.fq ref=$ref" >> Ribosomal_${file}.pbs ; done

# Paired, Corrected and Trimmed reads:
cut -f1 ListOfFiles | while read file; do f1="${MainDirectory}/Trimmomatic/CorTrimmed_${file}_R1.fq.gz"; f2="${MainDirectory}/Trimmomatic/CorTrimmed_${file}_R2.fq.gz"; echo "$toolsdir/bbmap/bbduk.sh t=$cpuless -Xmx${memless}g in=${f1} in2=${f2} outm=${MainDirectory}/RiboDepleted/riboRNA_${file}_R1.cor.fq.gz outm2=${MainDirectory}/RiboDepleted/riboRNA_${file}_R2.cor.fq.gz out=${MainDirectory}/RiboDepleted/CorTrimmed_${file}_R1.cor.fq out2=${MainDirectory}/RiboDepleted/CorTrimmed_${file}_R2.cor.fq ref=$ref" >> Ribosomal_${file}.pbs ; done

############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "===> Ribosomal reads depletion" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- command lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n2 Ribosomal_*.pbs >> $MainDirectory/Pipeline_Log

############################################
##       SUBMIT RIBODEPLETION JOBS        ##
############################################
# Populate Log file (header):
echo "Job Submission IDs:" >> ${MainDirectory}/PBSsub/SubJobs_Ribosomal

# Submit jobs and populate Log file (submission ID numbers):
ls Ribosomal_*.pbs | while read line; do echo "$line" >> ${MainDirectory}/PBSsub/SubJobs_Ribosomal ; qsub "$line"; done >> ${MainDirectory}/PBSsub/SubJobs_Ribosomal

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user" | grep Ribosomal
qstat | grep -c Ribosomal

############################################
##                   END                  ##
############################################

