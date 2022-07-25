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

step="READ TRIMMING"

echo "This script will run STEP 4 of the DENTAL pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Trim reads using Trimmomatic"

#  PIPELINE:
#  --------------------------------------------------------------
#   1 FASTQC PRE TRIM (quality-control before trimming)
#   2 RCORRECTOR (read correction)
#   3 FURpe (filtering uncorrectable reads)
#-->4 TRIMMOMATIC (trimming reads)
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

#   4 TRIMMOMATIC (trimming reads)
############################################
##             SOURCE OF TOOLS            ##
############################################

# Trimmomatic was obtained from:
# http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-Src-0.36.zip

# BBDuk is part of BBTools:
https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/

############################################
##           SOURCE OF ADAPTERS           ##
############################################

# WARNING!
# In order to define which adapter to trim (especially when sequencing details are not fully known) one can consult FastQC results and check for adapter contamination.
# These are examples of common adapters from http://docs.blast2go.com/user-manual/tools-(pro-feature)/fastq-quality-check/#FASTQQualityCheck-AdapterContent
# Illumina Universal Adapter: AGATCGGAAGAG (12bp)
# Nextera Transposase Sequence: CTGTCTCTTATA (12bp)
# Adapter sequences are available in the QIMRB cluster: /software/trimmomatic/trimmomatic-0.36/adapters/

############################################
##             TOOLS DIRECTORY            ##
############################################
# User-installed tools
# Local directory where all tools below are installed:
toolsdir="/working/lab_julietF/mainaB/Tools"

############################################
##            REQUIRED MODULES            ##
############################################
trimmomatic=trimmomatic/0.36
bbduk="${toolsdir}/bbmap/bbduk.sh"

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

# Create output directories:
mkdir ${MainDirectory}/Trimmomatic
mkdir ${MainDirectory}/BBduk


############################################
##            WRITE TO LOG FILE           ##
############################################

echo -e "\n\n" >> ${MainDirectory}/Pipeline_Log
echo "*******************************" >> ${MainDirectory}/Pipeline_Log
echo "Trimmomatic directory" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/Trimmomatic" >> ${MainDirectory}/Pipeline_Log
echo "What: Results of trimming reads with Trimmomatic" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) to reflect your requirements:

# Memory required (in GB):
memreq=25

# CPUs required (integer):
cpureq=24


############################################
##    CREATE PBS FILES FOR R-CORRECTOR    ##
############################################
# PBS files for this step are named Trimmomatic_FILE.pbs
# Cluster job names for this step are named Trimmomatic_FILE

# List all fastq files and create correspondent PBS files:
cut -f1 ListOfFiles | while read file; do cp ${MainDirectory}/Header.pbs Trimmomatic_${file}.pbs; done

# Populate PBS files (name of job):
cut -f1 ListOfFiles | while read file; do echo "#PBS -N Trimmomatic_${file}" >> Trimmomatic_${file}.pbs; done

# Populate PBS files (number of threads):
cut -f1 ListOfFiles | while read file; do echo "#PBS -r n" >> Trimmomatic_${file}.pbs; done

# Populate PBS files (time and resources for running the job):
cut -f1 ListOfFiles | while read file; do echo "#PBS -l mem=${memreq}GB,walltime=48:00:00,ncpus=${cpureq}" >> Trimmomatic_${file}.pbs; done

# Populate PBS files (send message when job is finished or aborted):
cut -f1 ListOfFiles | while read file; do echo "#PBS -m ae" >> Trimmomatic_${file}.pbs; done

# Populate PBS files (e-mail address for correspondence):
cut -f1 ListOfFiles | while read file; do echo "#PBS -M $email" >> Trimmomatic_${file}.pbs; done

# Populate PBS files (change directory):
cut -f1 ListOfFiles | while read file; do echo "cd $MainDirectory" >> Trimmomatic_${file}.pbs; done

# Populate PBS files (load module(s)):
cut -f1  ListOfFiles | while read file; do echo "module load $trimmomatic" >> Trimmomatic_${file}.pbs; done
cut -f1  ListOfFiles | while read file; do echo "module load $bbduk" >> Trimmomatic_${file}.pbs; done


# Populate PBS files (empty line):
cut -f1 ListOfFiles | while read file; do echo "" >> Trimmomatic_${file}.pbs; done


############################################
##       COMMAND LINE FOR TRIMMOMATIC     ##
############################################
# Populate PBS files (actual command line):

# WARNING!
# These two command lines will run Trimmomatic both for corrected and uncorrected reads, separately.
# If one wish to trim only corrected reads, comment or delete the first of the two command lines below.
# BBDuk is used here as a tool to measure overall read quality and produce graphs that can help find problematic elements (e.g. failed tiles, cycles, contaminations, etc).

# Trimming of uncorrected reads:
cut -f1 ListOfFiles | while read file; do f1="${MainDirectory}/${file}_R1.fq.gz"; f2="${MainDirectory}/${file}_R2.fq.gz"; echo "trimmomatic PE -phred33 -threads $cpureq ${f1} ${f2} ${MainDirectory}/Trimmomatic/Trimmed_${file}_R1.fq.gz ${MainDirectory}/Trimmomatic/Unpaired_Trimmed_${file}_R1.fq.gz ${MainDirectory}/Trimmomatic/Trimmed_${file}_R2.fq.gz ${MainDirectory}/Trimmomatic/Unpaired_Trimmed_${file}_R2.fq.gz ILLUMINACLIP:/software/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:75 HEADCROP:12" >> Trimmomatic_${file}.pbs ; done

# Trimming of corrected reads:
cut -f1 ListOfFiles | while read file; do f1="${MainDirectory}/FURpe/unfixrm_${file}_R1.cor.fq"; f2="${MainDirectory}/FURpe/unfixrm_${file}_R2.cor.fq"; echo "trimmomatic PE -phred33 -threads $cpureq ${f1} ${f2} ${MainDirectory}/Trimmomatic/CorTrimmed_${file}_R1.fq.gz ${MainDirectory}/Trimmomatic/Unpaired_CorTrimmed_${file}_R1.fq.gz ${MainDirectory}/Trimmomatic/CorTrimmed_${file}_R2.fq.gz ${MainDirectory}/Trimmomatic/Unpaired_CorTrimmed_${file}_R2.fq.gz ILLUMINACLIP:/software/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:75 HEADCROP:12" >> Trimmomatic_${file}.pbs ; done

# BBduk of uncorrected reads:
cut -f1 ListOfFiles | while read file; do f1="${MainDirectory}/${file}_R1.fq.gz"; f2="${MainDirectory}/${file}_R2.fq.gz"; echo "$bbduk t=$cpureq -Xmx${memreq}g in=${f1} in2=${f2} qtrim=rl trimq=30 minavgquality=20 minbasequality=10 minlen=75 outm=${MainDirectory}/BBduk/Failed_${file}_R1.fq.gz outm2=${MainDirectory}/BBduk/Failed_${file}_R2.fq.gz bhist=${MainDirectory}/BBduk/${file}_bhist.txt qhist=${MainDirectory}/BBduk/${file}_qhist.txt gchist=${MainDirectory}/BBduk/${file}_gchist.txt aqhist=${MainDirectory}/BBduk/${file}_aqhist.txt lhist=${MainDirectory}/BBduk/${file}_lhist.txt gcbins=auto" >> Trimmomatic_${file}.pbs ; done

# BBduk of corrected reads:
cut -f1 ListOfFiles | while read file; do f1="${MainDirectory}/Trimmomatic/CorTrimmed_${file}_R1.fq.gz"; f2="${MainDirectory}/Trimmomatic/CorTrimmed_${file}_R2.fq.gz"; echo "$bbduk t=$cpureq -Xmx${memreq}g in=${f1} in2=${f2} qtrim=rl trimq=30 minavgquality=20 minbasequality=10 minlen=75 outm=${MainDirectory}/BBduk/CorTrimmed_${file}_R1.cor.fq.gz outm2=${MainDirectory}/BBduk/CorTrimmed_${file}_R2.cor.fq.gz bhist=${MainDirectory}/BBduk/CorTrimmed_${file}_bhist.txt qhist=${MainDirectory}/BBduk/CorTrimmed_${file}_qhist.txt gchist=${MainDirectory}/BBduk/CorTrimmed_${file}_gchist.txt aqhist=${MainDirectory}/BBduk/CorTrimmed_${file}_aqhist.txt lhist=${MainDirectory}/BBduk/CorTrimmed_${file}_lhist.txt gcbins=auto" >> Trimmomatic_${file}.pbs ; done


############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (empty line):
echo "" >> ${MainDirectory}/Pipeline_Log

# Populate Log file (header):
echo "** Trimmomatic Command Lines **" >> ${MainDirectory}/Pipeline_Log

# Populate Log file (command lines used):
tail -n1 Trimmomatic_*.pbs >> ${MainDirectory}/Pipeline_Log

# Populate Log file (empty line):
echo "" >> ${MainDirectory}/Pipeline_Log

############################################
##        SUBMIT TRIMMOMATIC JOBS         ##
############################################
# Populate Log file (header):
echo "Trimmomatic Job Submission IDs:" >> ${MainDirectory}/PBSsub/SubJobs_Trimmomatic

# Submit jobs and populate Log file (submission ID numbers):
ls Trimmomatic_*.pbs | while read line; do echo "$line" >> ${MainDirectory}/PBSsub/SubJobs_Trimmomatic ; qsub "$line"; done >> ${MainDirectory}/PBSsub/SubJobs_Trimmomatic

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user" | grep Trimmomatic
qstat | grep -c Trimmomatic

############################################
##                   END                  ##
############################################


