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

step="MAPPING READS TO GENOME"

echo "This script will run STEP 8 of the DENTAL pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Mapping reads to human genome using Bowtie2"

#  PIPELINE:
#  --------------------------------------------------------------
#   1 FASTQC PRE TRIM (quality-control before trimming)
#   2 RCORRECTOR (read correction)
#   3 FURpe (filtering uncorrectable reads)
#   4 TRIMMOMATIC (trimming reads)
#   5 FASTQC POST TRIM (quality-control after trimming)
#   6 FILTERING RIBOSOMAL READS (filter rRNA with BBduk)
#   7 FIXING READ HEADERS (in-house script to fix formating)
#-->8 MAP TO GENOME (using Bowtie to flag contamination)
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

#   8 MAP TO GENOME (using Bowtie to flag contamination)
############################################
##            REQUIRED MODULES            ##
############################################
bowtie2=bowtie2/2.2.9
samtools=samtools/1.9

############################################
##           SET REFERENCE VARs           ##
############################################
# Main directory serving as repository for references:
referencedir=/working/lab_julietF/mainaB/ReferenceGenomes

# Bowtie2 index:
index=${referencedir}/GenomeGRCh38p7_bowtie2

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
mkdir $MainDirectory/GenAli

############################################
##            WRITE TO LOG FILE           ##
############################################

echo "Bowtie2 Genome Alignment directory" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/GenAli" >> ${MainDirectory}/Pipeline_Log
echo "What: Results for reads-to-genome alignment using Bowtie2" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) if needed:

# Memory required (in GB):
memreq=12

# CPUs required (integer):
cpureq=8

############################################
##      CREATE PBS FILES FOR BOWTIE2      ##
############################################

# PBS files for this step are named Bowtie_FILE.pbs
# Cluster job names for this step are named Bowtie_FILE

# List all fastq files and create correspondent PBS files:
cut -f1  ListOfFiles | while read file; do cp ${MainDirectory}/Header.pbs Bowtie_${file}.pbs; done

# Populate PBS files (name of job):
cut -f1  ListOfFiles | while read file; do echo "#PBS -N Bowtie_${file}" >> Bowtie_${file}.pbs; done

# Populate PBS files (number of threads):
cut -f1  ListOfFiles | while read file; do echo "#PBS -r n" >> Bowtie_${file}.pbs; done

# Populate PBS files (time and resources for running the job):
cut -f1  ListOfFiles | while read file; do echo "#PBS -l mem=${memreq}GB,walltime=08:00:00,ncpus=${cpureq}" >> Bowtie_${file}.pbs; done

# Populate PBS files (send message when job is finished or aborted):
cut -f1  ListOfFiles | while read file; do echo "#PBS -m ae" >> Bowtie_${file}.pbs; done

# Populate PBS files (e-mail address for correspondence):
cut -f1  ListOfFiles | while read file; do echo "#PBS -M $email" >> Bowtie_${file}.pbs; done

# Populate PBS files (empty line):
cut -f1  ListOfFiles | while read file; do echo "" >> Bowtie_${file}.pbs; done

# Populate PBS files (load module(s)):
cut -f1  ListOfFiles | while read file; do echo "module load $bowtie2" >> Bowtie_${file}.pbs; done

# Populate PBS files (load module(s)):
cut -f1  ListOfFiles | while read file; do echo "module load $samtools" >> Bowtie_${file}.pbs; done

# Populate PBS files (empty line):
cut -f1  ListOfFiles | while read file; do echo "" >> Bowtie_${file}.pbs; done

############################################
##        COMMAND LINE FOR BOWTIE2        ##
############################################
# Populate PBS files (actual command lines):

# Clean, corrected and formated reads:
cut -f1 ListOfFiles | while read file; do f1="${MainDirectory}/ReadyReads/CleanCorrected_${file}_R1.cor.fq.gz"; f2="${MainDirectory}/ReadyReads/CleanCorrected_${file}_R2.cor.fq.gz"; echo "bowtie2 --fast --no-unal -p${cpureq} -x ${index} -1 ${f1} -2 ${f2} | samtools view -@${cpureq} -bS - >> ${MainDirectory}/GenAli/CleanCorrected_${file}.bam" >> Bowtie_${file}.pbs ; done

cut -f1 ListOfFiles | while read file; do echo "samtools sort -@${cpureq} ${MainDirectory}/GenAli/CleanCorrected_${file}.bam -o ${MainDirectory}/GenAli/CleanCorrected_${file}.coordSorted.bam" >> Bowtie_${file}.pbs; done

cut -f1 ListOfFiles | while read file; do echo "samtools index -@${cpureq} ${MainDirectory}/GenAli/CleanCorrected_${file}.coordSorted.bam ${MainDirectory}/GenAli/CleanCorrected_${file}.coordSorted.bam.bai" >> Bowtie_${file}.pbs; done


############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "===> Genome alignment using Bowtie2" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- command lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n3 Bowtie*.pbs >> $MainDirectory/Pipeline_Log

############################################
##          SUBMIT BOWTIE2 JOBS           ##
############################################
# Populate Log file (header):
echo "Job Submission IDs:" >> ${MainDirectory}/PBSsub/SubJobs_Bowtie

# Submit jobs and populate Log file (submission ID numbers):
ls Bowtie_*.pbs | while read line; do echo "$line" >> ${MainDirectory}/PBSsub/SubJobs_Bowtie ; qsub "$line"; done >> ${MainDirectory}/PBSsub/SubJobs_Bowtie

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user" | grep Bowtie
qstat | grep -c Bowtie

############################################
##                   END                  ##
############################################

