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

step="TRANSCRIPTOME READ REPRESENTATION"

echo "This script will run STEP 14 of the DENTAL pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Align reads to assembly with Bowtie!"

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
#-->14 TRANSCRIPTOME READ REPRESENTATION
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

#   14 TRANSCRIPTOME READ REPRESENTATION
############################################
##              MAIN SOURCES              ##
############################################

# External sources:
# https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
# https://github.com/trinityrnaseq/trinityrnaseq/wiki/RNA-Seq-Read-Representation-by-Trinity-Assembly

############################################
##               INFORMATION              ##
############################################

###### WARNING! ######
#  User may omit the #
# last step. -Check! #
###### WARNING! ######

# If fragment size is known and user would like to skip the correspondent command line (bamPEFragmentSize), simply delete or comment the line and inform mean insert size below:
# fragmentsize=enternumberhere

############################################
##            REQUIRED MODULES            ##
############################################

bowtie2=bowtie2/2.2.9
samtools=samtools/1.9
python3=python/3.6.1
R=R/3.3.1

# Picard Tools:
picardtool="java -jar /software/picard/picard-tools-2.19.0/picard.jar"

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

echo "Alignment of reads to assembly for transcripts support" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/Bowtie" >> ${MainDirectory}/Pipeline_Log
echo "What: Assembly to reads alignments with Bowtie2" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) if needed:

# Memory required (in GB):
memreq=36

# CPUs required (integer):
cpureq=12
cpuless=10

############################################
##      CREATE PBS FILES FOR BOWTIE2      ##
############################################

# PBS files for this step are named Bowtie3nity_id_PROJECT.pbs
# Cluster job names for this step are named Bowtie3nity_id_PROJECT

# List all fastq files and create correspondent PBS files:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do cp ${MainDirectory}/Header.pbs Bowtie3nity_${id}_${project}.pbs; done; done

# Populate PBS files (name of job):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -N Bowtie3nity_${id}_${project}" >> Bowtie3nity_${id}_${project}.pbs; done; done

# Populate PBS files (number of threads):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -r n" >> Bowtie3nity_${id}_${project}.pbs; done; done

# Populate PBS files (time and resources for running the job):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -l mem=${memreq}GB,walltime=48:00:00,ncpus=${cpureq}" >> Bowtie3nity_${id}_${project}.pbs; done; done

# Populate PBS files (send message when job is finished or aborted):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -m ae" >> Bowtie3nity_${id}_${project}.pbs; done; done

# Populate PBS files (e-mail address for correspondence):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -M $email" >> Bowtie3nity_${id}_${project}.pbs; done; done

# Populate PBS files (empty line):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Bowtie3nity_${id}_${project}.pbs; done; done

# Populate PBS files (load module(s)):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $bowtie2" >> Bowtie3nity_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $samtools" >> Bowtie3nity_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $python3" >> Bowtie3nity_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $R" >> Bowtie3nity_${id}_${project}.pbs; done; done

# Populate PBS files (empty line):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Bowtie3nity_${id}_${project}.pbs; done; done


############################################
##        COMMAND LINE FOR BOWTIE2        ##
############################################
# Populate PBS files (actual command lines):

# Create new directory to store outputs:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "mkdir ${MainDirectory}/Bowtie_${id}_${project}" >> Bowtie3nity_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Bowtie3nity_${id}_${project}.pbs; done; done

# Build index from assemblies:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "bowtie2-build --threads ${cpuless} ${MainDirectory}/Trinity/Trinity_${id}_${project}.Trinity.fasta ${MainDirectory}/Bowtie_${id}_${project}/IndexTrinity_${id}_${project}" >> Bowtie3nity_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Bowtie3nity_${id}_${project}.pbs; done; done

# Run alignment with Bowtie:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do leftlist="${MainDirectory}/Trinity_${id}_${project}/insilico_read_normalization/left.norm.fq"; rightlist="${MainDirectory}/Trinity_${id}_${project}/insilico_read_normalization/right.norm.fq"; echo "bowtie2 -p ${cpuless} -q --no-unal -k 20 -x ${MainDirectory}/Bowtie_${id}_${project}/IndexTrinity_${id}_${project} -1 ${leftlist} -2 ${rightlist} | samtools view -@${cpuless} -Sb -o ${MainDirectory}/Bowtie_${id}_${project}/Bowtie_${id}_${project}.bam" >> Bowtie3nity_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Bowtie3nity_${id}_${project}.pbs; done; done

# Sort alignment with Samtools:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "samtools sort -@${cpuless} ${MainDirectory}/Bowtie_${id}_${project}/Bowtie_${id}_${project}.bam -o ${MainDirectory}/Bowtie_${id}_${project}/Bowtie_${id}_${project}.coordSorted.bam" >> Bowtie3nity_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Bowtie3nity_${id}_${project}.pbs; done; done

# Index alignment with Samtools:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "samtools index -@${cpuless} ${MainDirectory}/Bowtie_${id}_${project}/Bowtie_${id}_${project}.coordSorted.bam" >> Bowtie3nity_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Bowtie3nity_${id}_${project}.pbs; done; done

# Index assemblies with Samtools:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "samtools faidx ${MainDirectory}/Trinity/Trinity_${id}_${project}.Trinity.fasta" >> Bowtie3nity_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Bowtie3nity_${id}_${project}.pbs; done; done

# Calculate fragment sizes with PicardTools (faster):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "$picardtool CollectInsertSizeMetrics I=${MainDirectory}/Bowtie_${id}_${project}/Bowtie_${id}_${project}.coordSorted.bam O=${MainDirectory}/Bowtie_${id}_${project}/InsertSize_${id}_${project}.txt H=${MainDirectory}/Bowtie_${id}_${project}/InsertSizeHistograp_${id}_${project}.pdf INCLUDE_DUPLICATES=FALSE" >> Bowtie3nity_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Bowtie3nity_${id}_${project}.pbs; done; done


############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "===> Transcriptome Read Support" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- command lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n15 Bowtie3nity_*.pbs >> $MainDirectory/Pipeline_Log

############################################
##          SUBMIT BOWTIE2 JOBS           ##
############################################
# Populate Log file (header):
echo "Job Submission IDs:" >> ${MainDirectory}/PBSsub/SubJobs_Bowtie2

# Submit jobs and populate Log file (submission ID numbers):
ls Bowtie3nity_*.pbs | while read line; do echo "$line" >> ${MainDirectory}/PBSsub/SubJobs_Bowtie2 ; qsub "$line"; done >> ${MainDirectory}/PBSsub/SubJobs_Bowtie2

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user" | grep Bowtie3nity
qstat | grep -c Bowtie3nity

############################################
##                   END                  ##
############################################


