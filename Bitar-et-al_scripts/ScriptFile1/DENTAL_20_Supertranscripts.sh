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

step="MERGE TRANSCRIPTS TO SUPER TRANSCRIPTS"

echo "This script will run STEP 20 of the DENTAL pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Congregate transcripts and assign to 'supertranscrpts' (i.e. genes)"

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
#   18 ALIGN TRANSCRIPTS TO GENOME AND ASSESS SPLICING
#   19 FILTER OUT LOW SUPPORT TRANSCRIPTS
#-->20 MERGE TRANSCRIPTS TO SUPER TRANSCRIPTS
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

#   20 MERGE TRANSCRIPTS TO SUPER TRANSCRIPTS
############################################
##              MAIN SOURCES              ##
############################################


############################################
##               INFORMATION              ##
############################################


############################################
##            REQUIRED MODULES            ##
############################################

# Local directory where all tools below are installed:
toolsdir="/working/lab_julietF/mainaB/Tools"

# Trinity accessory scripts.
trinscripts="/software/trinityrnaseq/trinityrnaseq-2.8.4"
supertranscripts="${trinscripts}/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py"

# Modules:
python3=python/3.6.1
trinity=trinityrnaseq/2.8.4

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
mkdir ${MainDirectory}/SuperTranscripts

############################################
##            WRITE TO LOG FILE           ##
############################################

echo "FPKM support" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/Filter" >> ${MainDirectory}/Pipeline_Log
echo "What: Results of splicing determination" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) to reflect your requirements:

# Memory required (in GB):
memreq=32

# CPUs required (integer):
cpureq=16

############################################
##  CREATE PBS FILES FOR SUPERTRANSCRIPT  ##
############################################

# PBS files for this step are named Supertranscripts_id_PROJECT.pbs
# Cluster job names for this step are named Supertranscripts_id_PROJECT

# List all fastq files and create correspondent PBS files:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do cp ${MainDirectory}/Header.pbs Supertranscripts_${id}_${project}.pbs; done; done

# Populate PBS files (name of job):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -N Supertranscripts_${id}_${project}" >> Supertranscripts_${id}_${project}.pbs; done; done

# Populate PBS files (number of threads):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -r n" >> Supertranscripts_${id}_${project}.pbs; done; done

# Populate PBS files (time and resources for running the job):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -l mem=${memreq}GB,walltime=320:00:00,ncpus=${cpureq}" >> Supertranscripts_${id}_${project}.pbs; done; done

# Populate PBS files (send message when job is finished or aborted):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -m ae" >> Supertranscripts_${id}_${project}.pbs; done; done

# Populate PBS files (e-mail address for correspondence):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -M $email" >> Supertranscripts_${id}_${project}.pbs; done; done

# Populate PBS files (empty line):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Supertranscripts_${id}_${project}.pbs; done; done

# Populate PBS files (load module(s)):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $trinity" >> Supertranscripts_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $python3" >> Supertranscripts_${id}_${project}.pbs; done; done

# Populate PBS files (empty line):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Supertranscripts_${id}_${project}.pbs; done; done


############################################
##    COMMAND LINE FOR SUPERTRANSCRIPT    ##
############################################
# Populate PBS files (actual command line):

### HEADER file (create only once!)
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "grep '^>' ${MainDirectory}/Trinity/Trinity_${id}_${project}.Trinity.fasta >> ${MainDirectory}/SuperTranscripts/TrID_${id}_${project}" >> Supertranscripts_${id}_${project}.pbs; done; done


### RAW files

# Enter directory:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd ${MainDirectory}/SuperTranscripts" >> Supertranscripts_${id}_${project}.pbs; done; done

# Higher FPKM filter (1 // 5):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cat ${MainDirectory}/FPKMcutRaw_${id}_${project}/FPKM_${id}_${project}.spliced.fasta ${MainDirectory}/FPKMcutRaw_${id}_${project}/FPKM_${id}_${project}.unspliced.fasta >> ${MainDirectory}/SuperTranscripts/rawFPKM_${id}_${project}.fasta" >> Supertranscripts_${id}_${project}.pbs; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "$supertranscripts --trinity_fasta ${MainDirectory}/SuperTranscripts/rawFPKM_${id}_${project}.fasta --out_prefix rawSuper_${id}_${project}" >> Supertranscripts_${id}_${project}.pbs; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cat  ${MainDirectory}/SuperTranscripts/rawSuper_${id}_${project}.gtf"' | cut -f1 | grep "^TRINITY" | sort | uniq -c | cut -c1-8 | sort | uniq -c |'" sed '"'s/ \+ /\t/g'"' |  sed '"'s/^\t//g'"' >>  ${MainDirectory}/SuperTranscripts/rawSuper_${id}_${project}.graph" >> Supertranscripts_${id}_${project}.pbs; done; done


# Lower FPKM filter (0.5 // 3):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cat ${MainDirectory}/FPKMcutRaw_${id}_${project}/LowFPKM_${id}_${project}.spliced.fasta ${MainDirectory}/FPKMcutRaw_${id}_${project}/LowFPKM_${id}_${project}.unspliced.fasta >> ${MainDirectory}/SuperTranscripts/rawLowFPKM_${id}_${project}.fasta" >> Supertranscripts_${id}_${project}.pbs; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "$supertranscripts --trinity_fasta ${MainDirectory}/SuperTranscripts/rawLowFPKM_${id}_${project}.fasta --out_prefix rawLowSuper_${id}_${project}" >> Supertranscripts_${id}_${project}.pbs; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cat  ${MainDirectory}/SuperTranscripts/rawLowSuper_${id}_${project}.gtf"' | cut -f1 | grep "^TRINITY" | sort | uniq -c | cut -c1-8 | sort | uniq -c |'" sed '"'s/ \+ /\t/g'"' |  sed '"'s/^\t//g'"' >>  ${MainDirectory}/SuperTranscripts/rawLowSuper_${id}_${project}.graph" >> Supertranscripts_${id}_${project}.pbs; done; done


### FILTERED files

# Enter directory:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd ${MainDirectory}/SuperTranscripts" >> Supertranscripts_${id}_${project}.pbs; done; done

# Higher FPKM filter (1 // 5):

# Omit path information in fasta file:
#echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cat ${MainDirectory}/FPKMcutFiltered_${id}_${project}/FPKM_${id}_${project}.spliced.fasta ${MainDirectory}/FPKMcutFiltered_${id}_${project}/FPKM_${id}_${project}.unspliced.fasta >> ${MainDirectory}/SuperTranscripts/filteredFPKM_${id}_${project}.fasta" >> Supertranscripts_${id}_${project}.pbs; done; done

# With path information on fasta file:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cat ${MainDirectory}/FPKMcutFiltered_${id}_${project}/FPKM_${id}_${project}.spliced.fasta ${MainDirectory}/FPKMcutFiltered_${id}_${project}/FPKM_${id}_${project}.unspliced.fasta"' | grep "^>" | while read trid; do grep "${trid} len="'" ${MainDirectory}/SuperTranscripts/TrID_${id}_${project} >> ${MainDirectory}/SuperTranscripts/COL1_filteredFPKM_${id}_${project}; done; cat ${MainDirectory}/FPKMcutFiltered_${id}_${project}/FPKM_${id}_${project}.spliced.fasta ${MainDirectory}/FPKMcutFiltered_${id}_${project}/FPKM_${id}_${project}.unspliced.fasta"' | grep -v "^>"' ">> ${MainDirectory}/SuperTranscripts/COL2_filteredFPKM_${id}_${project}; paste ${MainDirectory}/SuperTranscripts/COL1_filteredFPKM_${id}_${project} ${MainDirectory}/SuperTranscripts/COL2_filteredFPKM_${id}_${project} | tr '\t' '\n' >> ${MainDirectory}/SuperTranscripts/filteredFPKM_${id}_${project}.fasta; wc -l ${MainDirectory}/SuperTranscripts/COL1_filteredFPKM_${id}_${project} ${MainDirectory}/SuperTranscripts/COL2_filteredFPKM_${id}_${project} >> ${MainDirectory}/SuperTranscripts/SizeLog_COLfiles; rm -rf ${MainDirectory}/SuperTranscripts/COL1_filteredFPKM_${id}_${project} ${MainDirectory}/SuperTranscripts/COL2_filteredFPKM_${id}_${project}" >> Supertranscripts_${id}_${project}.pbs; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "$supertranscripts --trinity_fasta ${MainDirectory}/SuperTranscripts/filteredFPKM_${id}_${project}.fasta --out_prefix filteredSuper_${id}_${project}" >> Supertranscripts_${id}_${project}.pbs; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cat  ${MainDirectory}/SuperTranscripts/filteredSuper_${id}_${project}.gtf"' | cut -f1 | grep "^TRINITY" | sort | uniq -c | cut -c1-8 | sort | uniq -c |'" sed '"'s/ \+ /\t/g'"' |  sed '"'s/^\t//g'"' >>  ${MainDirectory}/SuperTranscripts/filteredSuper_${id}_${project}.graph" >> Supertranscripts_${id}_${project}.pbs; done; done


# Lower FPKM filter (0.5 // 3):

# Omit path information in fasta file:
#echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cat ${MainDirectory}/FPKMcutFiltered_${id}_${project}/LowFPKM_${id}_${project}.spliced.fasta ${MainDirectory}/FPKMcutFiltered_${id}_${project}/LowFPKM_${id}_${project}.unspliced.fasta >> ${MainDirectory}/SuperTranscripts/filteredLowFPKM_${id}_${project}.fasta" >> Supertranscripts_${id}_${project}.pbs; done; done

# With path information on fasta file:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cat ${MainDirectory}/FPKMcutFiltered_${id}_${project}/LowFPKM_${id}_${project}.spliced.fasta ${MainDirectory}/FPKMcutFiltered_${id}_${project}/LowFPKM_${id}_${project}.unspliced.fasta"' | grep "^>" | while read trid; do grep "${trid} len="'" ${MainDirectory}/SuperTranscripts/TrID_${id}_${project} >> ${MainDirectory}/SuperTranscripts/COL1_filteredLowFPKM_${id}_${project}; done; cat ${MainDirectory}/FPKMcutFiltered_${id}_${project}/LowFPKM_${id}_${project}.spliced.fasta ${MainDirectory}/FPKMcutFiltered_${id}_${project}/LowFPKM_${id}_${project}.unspliced.fasta"' | grep -v "^>"' ">> ${MainDirectory}/SuperTranscripts/COL2_filteredLowFPKM_${id}_${project}; paste ${MainDirectory}/SuperTranscripts/COL1_filteredLowFPKM_${id}_${project} ${MainDirectory}/SuperTranscripts/COL2_filteredLowFPKM_${id}_${project} | tr '\t' '\n' >> ${MainDirectory}/SuperTranscripts/filteredLowFPKM_${id}_${project}.fasta; wc -l ${MainDirectory}/SuperTranscripts/COL1_filteredLowFPKM_${id}_${project} ${MainDirectory}/SuperTranscripts/COL2_filteredLowFPKM_${id}_${project} >> ${MainDirectory}/SuperTranscripts/SizeLog_COLfiles; rm -rf ${MainDirectory}/SuperTranscripts/COL1_filteredLowFPKM_${id}_${project} ${MainDirectory}/SuperTranscripts/COL2_filteredLowFPKM_${id}_${project}" >> Supertranscripts_${id}_${project}.pbs; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "$supertranscripts --trinity_fasta ${MainDirectory}/SuperTranscripts/filteredLowFPKM_${id}_${project}.fasta --out_prefix filteredLowSuper_${id}_${project}" >> Supertranscripts_${id}_${project}.pbs; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cat  ${MainDirectory}/SuperTranscripts/filteredLowSuper_${id}_${project}.gtf"' | cut -f1 | grep "^TRINITY" | sort | uniq -c | cut -c1-8 | sort | uniq -c |'" sed '"'s/ \+ /\t/g'"' |  sed '"'s/^\t//g'"' >>  ${MainDirectory}/SuperTranscripts/filteredLowSuper_${id}_${project}.graph" >> Supertranscripts_${id}_${project}.pbs; done; done



############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "===> Supertranscripts" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- command lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n20 Supertranscripts_*.pbs >> $MainDirectory/Pipeline_Log

############################################
##      SUBMIT SUPERTRANSCRIPT JOBS       ##
############################################
# Populate Log file (header):
echo "Job Submission IDs:" >> ${MainDirectory}/PBSsub/SubJobs_Supertranscripts

# Submit jobs and populate Log file (submission ID numbers):
ls Supertranscripts_*.pbs | while read file; do echo "$file" >> ${MainDirectory}/PBSsub/SubJobs_Supertranscripts ; qsub "$file"; done >> ${MainDirectory}/PBSsub/SubJobs_Supertranscripts

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user" | grep Supertranscripts
qstat | grep -c Supertranscripts

############################################
##                   END                  ##
############################################

