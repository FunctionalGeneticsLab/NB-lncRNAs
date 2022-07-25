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

step="GENERAL QUALITY STATISTICS"

echo "This script will run STEP 13 of the DENTAL pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Assess assembly!"

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
#-->13 GENERAL QUALITY STATISTICS
#   14 TRANSCRIPTOME READ REPRESENTATION
#   15 BUSCO (Benchmarking Universal Single-Copy Orthologs)
#   16 TransRate (Analyse de novo assembly)
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

#   13 GENERAL QUALITY STATISTICS
############################################
##              MAIN SOURCES              ##
############################################

# Trinity accessory scripts.
# [ to be used with "module load" ]

trinscripts="/software/trinityrnaseq/trinityrnaseq-2.8.4"

trinstats="${trinscripts}/util/TrinityStats.pl"
trinalign="${trinscripts}/util/align_and_estimate_abundance.pl"
trinfilter="${trinscripts}/util/filter_low_expr_transcripts.pl"
trinfilterold="/software/trinityrnaseq/trinityrnaseq-2.2.0/util/filter_fasta_by_rsem_values.pl"
trincdhit="$trinscripts/util/misc/filter_similar_seqs_expr_and_strand_aware.pl"
supertranscripts="${trinscripts}/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py"

############################################
##            REQUIRED MODULES            ##
############################################

trinity=trinityrnaseq/2.8.4
samtools=samtools/1.9

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
mkdir ${MainDirectory}/Trinity

############################################
##            WRITE TO LOG FILE           ##
############################################

echo "Trinity de novo assembly" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/Trinity" >> ${MainDirectory}/Pipeline_Log
echo "What: Results for de novo assembly using Trinity" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) if needed:

# Memory required (in GB):
memreq=6

# CPUs required (integer):
cpureq=1

############################################
##      CREATE PBS FILES FOR T-STATS      ##
############################################

# PBS files for this step are named Stats3nity_id_PROJECT.pbs
# Cluster job names for this step are named Stats3nity_id_PROJECT

# List all fastq files and create correspondent PBS files:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do cp ${MainDirectory}/Header.pbs Stats3nity_${id}_${project}.pbs; done; done

# Populate PBS files (name of job):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -N ${id}Stats3nity_${project}" >> Stats3nity_${id}_${project}.pbs; done; done

# Populate PBS files (number of threads):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -r n" >> Stats3nity_${id}_${project}.pbs; done; done

# Populate PBS files (time and resources for running the job):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -l mem=${memreq}GB,walltime=120:00:00,ncpus=${cpureq}" >> Stats3nity_${id}_${project}.pbs; done; done

# Populate PBS files (send message when job is finished or aborted):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -m ae" >> Stats3nity_${id}_${project}.pbs; done; done

# Populate PBS files (e-mail address for correspondence):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -M $email" >> Stats3nity_${id}_${project}.pbs; done; done

# Populate PBS files (empty line):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Stats3nity_${id}_${project}.pbs; done; done

# Populate PBS files (load module(s)):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $trinity" >> Stats3nity_${id}_${project}.pbs; done; done

# Populate PBS files (load module(s)):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $samtools" >> Stats3nity_${id}_${project}.pbs; done; done

# Populate PBS files (empty line):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Stats3nity_${id}_${project}.pbs; done; done


############################################
##        COMMAND LINE FOR T-STATS        ##
############################################
# Populate PBS files (actual command lines):

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "$trinstats ${MainDirectory}/Trinity/Trinity_${id}_${project}.Trinity.fasta >> ${MainDirectory}/Trinity/TrinityStats_${id}_${project}" >> Stats3nity_${id}_${project}.pbs; done; done


############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "===> Statistics of Raw Assembly" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- command lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n1 Stats3nity_*.pbs >> $MainDirectory/Pipeline_Log

############################################
##          SUBMIT TRINITY JOBS           ##
############################################
# Populate Log file (header):
echo "Job Submission IDs:" >> ${MainDirectory}/PBSsub/SubJobs_Stats3nity

# Submit jobs and populate Log file (submission ID numbers):
ls Stats3nity_*.pbs | while read line; do echo "$line" >> ${MainDirectory}/PBSsub/SubJobs_Stats3nity ; qsub "$line"; done >> ${MainDirectory}/PBSsub/SubJobs_Stats3nity

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user" | grep Stats3nity
qstat | grep -c Stats3nity

############################################
##                   END                  ##
############################################


