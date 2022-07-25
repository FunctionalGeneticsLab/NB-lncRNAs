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

step="TransRate (Analyse de novo assembly)"

echo "This script will run STEP 16 of the DENTAL pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Assess the quality of the assembly using Transrate!"

#  PIPELINE:
#  --------------------------------------------------------------
#   1 FASTQC PRE TRIM (quality-control before trimming)
#   2 RCORRECTOR (read correction)
#   3 FURpe (filtering uncorrectable reads)
#   4 TRIMMOMATIC (trimming reads)
#   5 FASTQC POST TRIM (quality-control after trimming)
#   6 FILTERING RIBOSOMAL READS (filter rRNA with BBduk)
#   7 FIXING READ HEADERS (in-house script to fix formating)
#   8 MAP TO GENOME (using Bowtie to flag contamination)
#   9 INFER EXPERIMENTAL DESIGN (check strandness using RSeQC)
#   10 FASTQC PRE TRINITY (final quality-control before Trinity)
#  --------------------------------------------------------------
#   11 ASSESS REPLICATES CORRELATION
#   12 TRINITY
#   13 GENERAL QUALITY STATISTICS
#   14 TRANSCRIPTOME READ REPRESENTATION
#   15 BUSCO (Benchmarking Universal Single-Copy Orthologs)
#-->16 TransRate (Analyse de novo assembly)
#   17 Detonate (DE novo TranscriptOme rNa-seq Assembly with or without the Truth Evaluation)
#   18 ALIGN TRANSCRIPTS TO GENOME AND ASSESS SPLICING
#   19 FILTER OUT LOW SUPPORT TRANSCRIPTS
#   20 MERGE TRANSCRIPTS TO SUPER TRANSCRIPTS
#  --------------------------------------------------------------
#   21 RETRIEVE ANNOTATION (to rule out batch effects)
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

#   16 TransRate (Analyse de novo assembly)
############################################
##              MAIN SOURCES              ##
############################################

# TransRate (pre-compiled version) was obtained from:
# https://bintray.com/artifact/download/blahah/generic/transrate-1.0.3-linux-x86_64.tar.gz
# /working/lab_julietF/mainaB/Tools/transrate-1.0.3-linux-x86_64/transrate --install-deps ref

############################################
##               INFORMATION              ##
############################################

# About fastq-pair table size calibration:
# The size of the hash table may be altered by the user, and a larger table, while consuming more memory and taking slightly longer to instantiate will reduce look-up time to determine whether an identifier is in the hash. The default table size (100,003) is designed to enhance even distribution of elements through the table, and provides linear complexity even for large sequence files.


############################################
##            REQUIRED MODULES            ##
############################################

# Local directory where all tools below are installed:
toolsdir="/working/lab_julietF/mainaB/Tools"

fastqpair="${toolsdir}/fastq-pair/build/fastq_pair"
#transrate="${toolsdir}/transrate-1.0.3-linux-x86_64/transrate"
transrate="${toolsdir}/transrate-1.0.1-linux-x86_64/transrate"

# Annotation files for human transcriptome.
# Genome build GRCh38.p2 updated in 2015.
# Main directory serving as repository for references:
referencedir=/working/lab_julietF/mainaB/ReferenceGenomes
transgtf=${referencedir}/TranscriptomeGRCh38rel79_ERCC.gtf
transplusgtf=${referencedir}/TranscriptomeGRCh38rel79_plus_mencRNAs.gtf
transcriptome=${referencedir}/TranscriptomeGRCh38p7_ERCC.fa

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
mkdir ${MainDirectory}/Transrate

############################################
##            WRITE TO LOG FILE           ##
############################################

echo "Transrate analysis" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/Transrate" >> ${MainDirectory}/Pipeline_Log
echo "What: Results for assembly quality assessment" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) if needed:

# Memory required (in GB):
memreq=80

# CPUs required (integer):
cpureq=32

############################################
##     CREATE PBS FILES FOR TRANSRATE     ##
############################################
# PBS files for this step are named TransRate_id_PROJECT.pbs
# Cluster job names for this step are named TransRate_id_PROJECT

# List all fastq files and create correspondent PBS files:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do cp ${MainDirectory}/Header.pbs TransRate_${id}_${project}.pbs; done; done

# Populate PBS files (name of job):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -N TransRate_${id}_${project}" >> TransRate_${id}_${project}.pbs; done; done

# Populate PBS files (number of threads):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -r n" >> TransRate_${id}_${project}.pbs; done; done

# Populate PBS files (time and resources for running the job):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -l mem=${memreq}GB,walltime=120:00:00,ncpus=${cpureq}" >> TransRate_${id}_${project}.pbs; done; done

# Populate PBS files (send message when job is finished or aborted):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -m ae" >> TransRate_${id}_${project}.pbs; done; done

# Populate PBS files (e-mail address for correspondence):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -M $email" >> TransRate_${id}_${project}.pbs; done; done

# Populate PBS files (empty line):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> TransRate_${id}_${project}.pbs; done; done

# Populate PBS files (empty line):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> TransRate_${id}_${project}.pbs; done; done


############################################
##       COMMAND LINE FOR TRANSRATE       ##
############################################
# Populate PBS files (actual command line):

echo ""; echo "Please wait while the program calculates the table sizes..."; echo ""

# Create output directories:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "mkdir ${MainDirectory}/Transrate_${id}_${project}" >> TransRate_${id}_${project}.pbs; done; done

# Properly pair all reads with fastq-pair:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do tablesize=`wc -l ${MainDirectory}/Trinity_${id}_${project}/insilico_read_normalization/left.norm.fq | cut -d' ' -f1`; echo "--> Table size:$tablesize"; echo "" ;echo "$fastqpair -t $tablesize ${MainDirectory}/Trinity_${id}_${project}/insilico_read_normalization/left.norm.fq ${MainDirectory}/Trinity_${id}_${project}/insilico_read_normalization/right.norm.fq" >> TransRate_${id}_${project}.pbs; done; done

# Run TransRate:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "$transrate --assembly ${MainDirectory}/Trinity/Trinity_${id}_${project}.Trinity.fasta --left ${MainDirectory}/Trinity_${id}_${project}/insilico_read_normalization/left.norm.fq.paired.fq --right ${MainDirectory}/Trinity_${id}_${project}/insilico_read_normalization/right.norm.fq.paired.fq --reference ${transcriptome} --threads ${cpureq} --output ${MainDirectory}/Transrate_${id}_${project}" >> TransRate_${id}_${project}.pbs; done; done


############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "===> Transrate assembly quality" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- command lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n3 TransRate_*.pbs >> $MainDirectory/Pipeline_Log

############################################
##         SUBMIT TRANSRATE  JOBS         ##
############################################
# Populate Log file (header):
echo "Job Submission IDs:" >> ${MainDirectory}/PBSsub/SubJobs_Transrate

# Submit jobs and populate Log file (submission ID numbers):
ls TransRate_*.pbs | while read line; do echo "$line" >> ${MainDirectory}/PBSsub/SubJobs_Transrate ; qsub "$line"; done >> ${MainDirectory}/PBSsub/SubJobs_Transrate

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user" | grep TransRate
qstat | grep -c TransRate

############################################
##                   END                  ##
############################################

