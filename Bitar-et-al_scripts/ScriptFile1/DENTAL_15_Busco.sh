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

step="BUSCO (Benchmarking Universal Single-Copy Orthologs)"

echo "This script will run STEP 15 of the DENTAL pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Assess the compleness of the assembly using BUSCO!"

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
#-->15 BUSCO (Benchmarking Universal Single-Copy Orthologs)
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

#   15 BUSCO (Benchmarking Universal Single-Copy Orthologs)
############################################
##              MAIN SOURCES              ##
############################################

## BUSCO
# BUSCO is available at the QIMRB cluster and Augustus was prepared using:
# tar -xvf /software/busco/busco-20161219/augustus-3.2.2.config.tgz

# BUSCO lineage data were obtained from:
# https://busco-archive.ezlab.org/v3/datasets


############################################
##               INFORMATION              ##
############################################

# WARNING!
# If needed (not previously done), untar augustus in the home folder:
# tar -xvf /software/busco/busco-20161219/augustus-3.2.2.config.tgz

# If needed (if not previously done), download lineage data from Busco (e.g. eukaryota, vertebrata or mammalia):
# https://busco-archive.ezlab.org/v3/datasets/eukaryota_odb9.tar.gz
# https://busco-archive.ezlab.org/v3/datasets/vertebrata_odb9.tar.gz
# https://busco-archive.ezlab.org/v3/datasets/mammalia_odb9.tar.gz

############################################
##            REQUIRED MODULES            ##
############################################

# Set species name:
species="human"

python3=python/3.6.1
busco=busco/20161119
oldblast=ncbi-blast/2.2.31+

# Local directory where all tools below are installed:
toolsdir="/working/lab_julietF/mainaB/Tools"
buscolineages="${toolsdir}/BuscoLineages"

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

echo "BUSCO analysis" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/BUSCO" >> ${MainDirectory}/Pipeline_Log
echo "What: Results for assembly completeness" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) if needed:

# Memory required (in GB):
memreq=24

# CPUs required (integer):
cpureq=1

############################################
##       CREATE PBS FILES FOR BUSCO       ##
############################################

# PBS files for this step are named Busco_id_PROJECT.pbs
# Cluster job names for this step are named Busco_id_PROJECT

# List all fastq files and create correspondent PBS files:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do cp ${MainDirectory}/Header.pbs Busco_${id}_${project}.pbs; done; done

# Populate PBS files (name of job):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -N Busco_${id}_${project}" >> Busco_${id}_${project}.pbs; done; done

# Populate PBS files (number of threads):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -r n" >> Busco_${id}_${project}.pbs; done; done

# Populate PBS files (time and resources for running the job):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -l mem=${memreq}GB,walltime=120:00:00,ncpus=${cpureq}" >> Busco_${id}_${project}.pbs; done; done

# Populate PBS files (send message when job is finished or aborted):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -m ae" >> Busco_${id}_${project}.pbs; done; done

# Populate PBS files (e-mail address for correspondence):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -M $email" >> Busco_${id}_${project}.pbs; done; done

# Populate PBS files (empty line):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Busco_${id}_${project}.pbs; done; done

# Populate PBS files (load module(s)):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $python3" >> Busco_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $busco" >> Busco_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $oldblast" >> Busco_${id}_${project}.pbs; done; done

# Populate PBS files (empty line):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Busco_${id}_${project}.pbs; done; done


############################################
##         COMMAND LINE FOR BUSCO         ##
############################################
# Populate PBS files (actual command lines):

# Enter directory:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd $MainDirectory" >> Busco_${id}_${project}.pbs; done; done

# Run BUSCO for eukaryota:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "BUSCO.py -i ${MainDirectory}/Trinity/Trinity_${id}_${project}.Trinity.fasta -l ${buscolineages}/eukaryota_odb9 -o BuscoEukaryota_${id}_${project} -m transcriptome -c ${cpureq} -t ${MainDirectory}/BuscoTempE_${id}_${project} -sp ${species}" >> Busco_${id}_${project}.pbs; done; done

# Run BUSCO for mammalia:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "BUSCO.py -i ${MainDirectory}/Trinity/Trinity_${id}_${project}.Trinity.fasta -l ${buscolineages}/mammalia_odb9 -o BuscoMammalia_${id}_${project} -m transcriptome -c ${cpureq} -t ${MainDirectory}/BuscoTempM_${id}_${project} -sp ${species}" >> Busco_${id}_${project}.pbs; done; done


############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "===> BUSCO assembly completeness" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- command lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n3 Busco_*.pbs >> $MainDirectory/Pipeline_Log

############################################
##           SUBMIT BUSCO  JOBS           ##
############################################
# Populate Log file (header):
echo "Job Submission IDs:" >> ${MainDirectory}/PBSsub/SubJobs_Busco

# Submit jobs and populate Log file (submission ID numbers):
ls Busco_*.pbs | while read line; do echo "$line" >> ${MainDirectory}/PBSsub/SubJobs_Busco ; qsub "$line"; done >> ${MainDirectory}/PBSsub/SubJobs_Busco

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user" | grep Busco
qstat | grep -c Busco

############################################
##                   END                  ##
############################################


