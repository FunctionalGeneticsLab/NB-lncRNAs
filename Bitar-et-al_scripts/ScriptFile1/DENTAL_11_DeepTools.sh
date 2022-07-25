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

step="ASSESS REPLICATES CORRELATION"

echo "This script will run STEP 11 of the DENTAL pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Rule out significant batch effects by assessing replicates correlation."

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
#-->11 ASSESS REPLICATES CORRELATION
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

#   11 ASSESS REPLICATES CORRELATION (to rule out batch effects)
############################################
##             SOURCE OF TOOLS            ##
############################################

# FastQC was obtained from:
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc

############################################
##            REQUIRED MODULES            ##
############################################
python3=python/3.6.1

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
mkdir ${MainDirectory}/DeepTools

############################################
##            WRITE TO LOG FILE           ##
############################################

echo "DeepTools replicate correlation" >> ${MainDirectory}/Pipeline_Log
echo "Where: ${MainDirectory}/DeepTools" >> ${MainDirectory}/Pipeline_Log
echo "What: Results for replicate correlation using DeepTools" >> ${MainDirectory}/Pipeline_Log
echo "" >> ${MainDirectory}/Pipeline_Log

############################################
##              SET RESOURCES             ##
############################################
# Edit the following line(s) if needed:

# Memory required (in GB):
memreq=26

# CPUs required (integer):
cpureq=8

############################################
##     CREATE PBS FILES FOR DEEPTOOLS     ##
############################################
# PBS files for this project are named deeptools_PROJECT.pbs
# Cluster job names for this project are named deeptools_PROJECT

project=`cut -d'_' -f1 ListOfFiles | sort | uniq`

# List all fastq files and create correspondent PBS files:
cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do cp ${MainDirectory}/Header.pbs deeptools_${project}.pbs; done

# Populate PBS files (name of job):
cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -N deeptools_${project}" >> deeptools_${project}.pbs; done

# Populate PBS files (number of threads):
cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -r n" >> deeptools_${project}.pbs; done

# Populate PBS files (time and resources for running the job):
cut -d'_' -f1 ListOfFiles | sort | uniq  | while read project; do echo "#PBS -l mem=${memreq}GB,walltime=08:00:00,ncpus=${cpureq}" >> deeptools_${project}.pbs; done

# Populate PBS files (send message when job is finished or aborted):
cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -m ae" >> deeptools_${project}.pbs; done

# Populate PBS files (e-mail address for correspondence):
cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -M $email" >> deeptools_${project}.pbs; done

# Populate PBS files (empty line):
cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> deeptools_${project}.pbs; done

# Populate PBS files (load module(s)):
cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $python3" >> deeptools_${project}.pbs; done

# Populate PBS files (empty line):
cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> deeptools_${project}.pbs; done


############################################
##       COMMAND LINE FOR DEEPTOOLS       ##
############################################
# Populate PBS files (actual command lines):

# List colors for DeepTools (in hex code or actual names):
colors="#61bcc6;#ade66c;#f11e71;#9cd5d8;#6213fe;#e17ebb;#f7dc02;#a14b74;#a66fc7;#de08e3;#be1f5f;#1838f9;#d6fd2e;#16a184;#5b1119;#b9f9f7;#57e148;#0c4b0e;#a88f99;#4b8623;#0c170d;#cfc154;#424cab;#669dd7;#3a06a8;#ec5a5e;#56dac1;#e4d2d6;#262f4f;#b3d17e"
colors="black;grey;lightgrey;palevioletred;firebrick;crimson;orangered;saddlebrown;bisque;navajowhite;orchid;magenta;darkkhaki;olive;darkolivegreen;darkseagreen;darkgreen;mediumseagreen;aquamarine;lightcyan;darkcyan;cadetblue;lightskyblue;lightslategrey;royalblue;darkblue;darkslateblue;indigo;plum;fuchsia;hotpink;lightpink"

# Generate list of colors for PCA plot:
echo ${colors} | tr ';' '\n' >> ${MainDirectory}/DeepTools/ColorPallete

# Generate color scheme for PCA plot:
project=`cut -d'_' -f1 ListOfFiles | sort | uniq`; labels=`grep -Pe "^${project}_" ListOfFiles | cut -d'_' -f2 | tr '\n' ' '`; echo ${labels} | tr ' ' '\n' | uniq >> ${MainDirectory}/DeepTools/tempLabels; samplenum=`cat ${MainDirectory}/DeepTools/tempLabels | uniq | wc -l`; head -n${samplenum} ${MainDirectory}/DeepTools/ColorPallete >> ${MainDirectory}/DeepTools/tempColors; paste ${MainDirectory}/DeepTools/tempLabels ${MainDirectory}/DeepTools/tempColors >> ${MainDirectory}/DeepTools/LabelColors_${project}; echo ${labels} | tr ' ' '\n' | while read sample; do grep -Pe "^${sample}\t" ${MainDirectory}/DeepTools/LabelColors_${project} | cut -f2 >> ${MainDirectory}/DeepTools/Colist_${project}; rm -rf ${MainDirectory}/DeepTools/temp*; done

# Use cleaned and Corrected paired-end reads previously aligned to genome.

# Generate BAM comparison summary:
project=`cut -d'_' -f1 ListOfFiles | sort | uniq`; labels=`grep -Pe "^${project}_" ListOfFiles | cut -d'_' -f2 | tr '\n' ' '`; ls ${MainDirectory}/GenAli/CleanCorrected_*.coordSorted.bam >> ${MainDirectory}/DeepTools/BamList_${project}; bams=`cat ${MainDirectory}/DeepTools/BamList_${project} | tr '\n' ' '`; echo "multiBamSummary bins --bamfiles $bams --labels $labels --outFileName ${MainDirectory}/DeepTools/Correlation_${project}.npz --maxFragmentLength 280 --numberOfProcessors ${cpureq}" >> deeptools_${project}.pbs

# Generate Spearman correlation:
project=`cut -d'_' -f1 ListOfFiles | sort | uniq`; labels=`grep -Pe "^${project}_" ListOfFiles | cut -d'_' -f2 | tr '\n' ' '`; echo "plotCorrelation --corData ${MainDirectory}/DeepTools/Correlation_${project}.npz --corMethod spearman --whatToPlot heatmap --labels $labels --plotFile ${MainDirectory}/DeepTools/SpearmanCorrelation_${project}.pdf --plotTitle SpearmanCorrelation_${project} --outFileCorMatrix ${MainDirectory}/DeepTools/SpearmanCorrelation_${project}.matrix --plotNumbers" >> deeptools_${project}.pbs

# Generate Pearson correlation:
project=`cut -d'_' -f1 ListOfFiles | sort | uniq`; labels=`grep -Pe "^${project}_" ListOfFiles | cut -d'_' -f2 | tr '\n' ' '`; echo "plotCorrelation --corData ${MainDirectory}/DeepTools/Correlation_${project}.npz --corMethod pearson --whatToPlot heatmap --labels $labels --plotFile ${MainDirectory}/DeepTools/PearsonCorrelation_${project}.pdf --plotTitle PearsonCorrelation_${project} --outFileCorMatrix ${MainDirectory}/DeepTools/PearsonCorrelation_${project}.matrix --plotNumbers" >> deeptools_${project}.pbs

# Generate PCA plot:
project=`cut -d'_' -f1 ListOfFiles | sort | uniq`; labels=`grep -Pe "^${project}_" ListOfFiles | cut -d'_' -f2 | tr '\n' ' '`; colors=`cat ${MainDirectory}/DeepTools/Colist_${project} | tr '\n' ' '`; echo "plotPCA --corData ${MainDirectory}/DeepTools/Correlation_${project}.npz --plotHeight 20 --labels $labels --colors ${colors} --plotFile ${MainDirectory}/DeepTools/PCAplot_${project}.pdf --plotTitle PCAplot_${project} --log2" >> deeptools_${project}.pbs


############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "===> Replicate correlation with DeepTools" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- command lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail -n4 deeptools*.pbs >> $MainDirectory/Pipeline_Log

############################################
##         SUBMIT DEEPTOOLS JOBS          ##
############################################
# Populate Log file (header):
echo "Job Submission IDs:" >> ${MainDirectory}/PBSsub/SubJobs_DeepTools

# Submit jobs and populate Log file (submission ID numbers):
ls deeptools*.pbs | while read line; do echo "$line" >> ${MainDirectory}/PBSsub/SubJobs_DeepTools ; qsub "$line"; done >> ${MainDirectory}/PBSsub/SubJobs_DeepTools

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user" | grep deeptools
qstat | grep -c deeptools

############################################
##                   END                  ##
############################################

