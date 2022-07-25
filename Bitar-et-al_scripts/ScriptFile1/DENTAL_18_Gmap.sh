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

step="ALIGN TRANSCRIPTS TO GENOME AND ASSESS SPLICING"

echo "This script will run STEP 18 of the DENTAL pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Align transcripts to genome and assess splicing!"

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
#   16 TransRate (Analyse de novo assembly)
#   17 Detonate (DE novo TranscriptOme rNa-seq Assembly with or without the Truth Evaluation)
#-->18 ALIGN TRANSCRIPTS TO GENOME AND ASSESS SPLICING
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

#   18 ALIGN TRANSCRIPTS TO GENOME AND ASSESS SPLICING
############################################
##              MAIN SOURCES              ##
############################################

# External sources:
# Information on human intron length distribution (selected cutoff was 50bp).
# https://www.researchgate.net/publication/7498905_An_analysis_on_gene_architecture_in_human_and_mouse_genomes
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0175393
# https://www.biostars.org/p/89581/

############################################
##               INFORMATION              ##
############################################

# Separate spliced x unspliced:
# The 'N' tag in BAM format represents skipped region from the reference. So if a read doesn't have a continuous alignment or a large reference region is skipped from the alignment then that portion of reference genome will be depicted in BAM file using 'N' tag in the sixth column. This cannot guarantee, however, that both aligned portions of the reads are mapped to exons. The aligned fragments can be covering exon-intron, intron-exon or intron-intron regions as well. But irrespective to RNA-seq reads, most probably these are exon-exon.

############################################
##            REQUIRED MODULES            ##
############################################

# Local directory where all tools below are installed:
toolsdir="/working/lab_julietF/mainaB/Tools"
gmap="${toolsdir}/gmap-2020-09-12/src/gmap"

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
mkdir ${MainDirectory}/Gmap

############################################
##            WRITE TO LOG FILE           ##
############################################

echo "Gmap alignment" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/Gmap" >> ${MainDirectory}/Pipeline_Log
echo "What: Results of splicing determination" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) if needed:

# Memory required (in GB):
memreq=32

# CPUs required (integer):
cpureq=16

############################################
##        CREATE PBS FILES FOR GMAP       ##
############################################
# PBS files for this step are named Gmap_id_PROJECT.pbs
# Cluster job names for this step are named Gmap_id_PROJECT

# List all fastq files and create correspondent PBS files:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do cp ${MainDirectory}/Header.pbs Gmap_${id}_${project}.pbs; done; done

# Populate PBS files (name of job):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -N Gmap_${id}_${project}" >> Gmap_${id}_${project}.pbs; done; done

# Populate PBS files (number of threads):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -r n" >> Gmap_${id}_${project}.pbs; done; done

# Populate PBS files (time and resources for running the job):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -l mem=${memreq}GB,walltime=120:00:00,ncpus=${cpureq}" >> Gmap_${id}_${project}.pbs; done; done

# Populate PBS files (send message when job is finished or aborted):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -m ae" >> Gmap_${id}_${project}.pbs; done; done

# Populate PBS files (e-mail address for correspondence):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -M $email" >> Gmap_${id}_${project}.pbs; done; done

# Populate PBS files (empty line):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Gmap_${id}_${project}.pbs; done; done

# Populate PBS files (empty line):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Gmap_${id}_${project}.pbs; done; done


############################################
##          COMMAND LINE FOR GMAP         ##
############################################
# Populate PBS files (actual command line):

# Divide Gmap in 2 (for raw and filtered assemblies):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do cat Gmap_${id}_${project}.pbs | sed 's/PBS -N Gmap_/PBS -N GmapFiltered_/g' >> GmapFiltered_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do cat Gmap_${id}_${project}.pbs | sed 's/PBS -N Gmap_/PBS -N GmapRaw_/g' >> GmapRaw_${id}_${project}.pbs; done; done

# Create output directories:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "mkdir ${MainDirectory}/GmapFiltered_${id}_${project}" >> GmapFiltered_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "mkdir ${MainDirectory}/GmapRaw_${id}_${project}" >> GmapRaw_${id}_${project}.pbs; done; done

# Run Gmap:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "$gmap -n 0 -t ${cpureq} -B 5 --dir=${referencedir} --db=GenomeGRCh38p7_gmap ${MainDirectory}/Transrate_${id}_${project}/Trinity_${id}_${project}.Trinity/good.Trinity_${id}_${project}.Trinity.fasta --min-intronlength=50 --max-intronlength-middle=200000 --min-identity=0.9 --format=samse > ${MainDirectory}/GmapFiltered_${id}_${project}/Transrate_${id}_${project}.genome.sam" >> GmapFiltered_${id}_${project}.pbs; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "$gmap -n 0 -t ${cpureq} -B 5 --dir=${referencedir} --db=GenomeGRCh38p7_gmap ${MainDirectory}/Trinity/Trinity_${id}_${project}.Trinity.fasta --min-intronlength=50 --max-intronlength-middle=200000 --min-identity=0.9 --format=samse > ${MainDirectory}/GmapRaw_${id}_${project}/Trinity_${id}_${project}.genome.sam" >> GmapRaw_${id}_${project}.pbs; done; done

# Divide transcripts in categories (spliced, unspliced, low-quality alignment and unaligned):

# Filtered Assembly:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cat ${MainDirectory}/GmapFiltered_${id}_${project}/Transrate_${id}_${project}.genome.sam | awk '"'$1 ~ /@/ || ($6 !~ /N/ && $6 !~ /*/)'"' >> ${MainDirectory}/GmapFiltered_${id}_${project}/Transrate_${id}_${project}.unspliced.sam" >> GmapFiltered_${id}_${project}.pbs; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cat ${MainDirectory}/GmapFiltered_${id}_${project}/Transrate_${id}_${project}.genome.sam | awk '"'$1 ~ /@/ || ($6 ~ /N/)'"' >> ${MainDirectory}/GmapFiltered_${id}_${project}/Transrate_${id}_${project}.spliced.sam" >> GmapFiltered_${id}_${project}.pbs; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cat ${MainDirectory}/GmapFiltered_${id}_${project}/Transrate_${id}_${project}.genome.sam | awk '"'$1 ~ /@/ || ($6 ~ /*/)'"' >> ${MainDirectory}/GmapFiltered_${id}_${project}/Transrate_${id}_${project}.unaligned.sam" >> GmapFiltered_${id}_${project}.pbs; done; done

# Raw Assembly:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cat ${MainDirectory}/GmapRaw_${id}_${project}/Trinity_${id}_${project}.genome.sam | awk '"'$1 ~ /@/ || ($6 !~ /N/ && $6 !~ /*/)'"' >> ${MainDirectory}/GmapRaw_${id}_${project}/Trinity_${id}_${project}.unspliced.sam" >> GmapRaw_${id}_${project}.pbs; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cat ${MainDirectory}/GmapRaw_${id}_${project}/Trinity_${id}_${project}.genome.sam | awk '"'$1 ~ /@/ || ($6 ~ /N/)'"' >> ${MainDirectory}/GmapRaw_${id}_${project}/Trinity_${id}_${project}.spliced.sam" >> GmapRaw_${id}_${project}.pbs; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cat ${MainDirectory}/GmapRaw_${id}_${project}/Trinity_${id}_${project}.genome.sam | awk '"'$1 ~ /@/ || ($6 ~ /*/)'"' >> ${MainDirectory}/GmapRaw_${id}_${project}/Trinity_${id}_${project}.unaligned.sam" >> GmapRaw_${id}_${project}.pbs; done; done


############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "===> Gmap splicing assessment" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- command lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n5 Gmap*.pbs >> $MainDirectory/Pipeline_Log

############################################
##            SUBMIT GMAP JOBS            ##
############################################
# Populate Log file (header):
echo "Job Submission IDs:" >> ${MainDirectory}/PBSsub/SubJobs_Gmap

# Submit jobs and populate Log file (submission ID numbers):
ls GmapFiltered_*.pbs | while read file; do echo "$file" >> ${MainDirectory}/PBSsub/SubJobs_Gmap ; qsub "$file"; done >> ${MainDirectory}/PBSsub/SubJobs_Gmap
ls GmapRaw_*.pbs | while read file; do echo "$file" >> ${MainDirectory}/PBSsub/SubJobs_Gmap ; qsub "$file"; done >> ${MainDirectory}/PBSsub/SubJobs_Gmap

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user" | grep Gmap
qstat | grep -c Gmap

############################################
##                   END                  ##
############################################

