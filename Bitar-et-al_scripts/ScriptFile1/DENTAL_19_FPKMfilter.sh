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

step="FILTER OUT LOW SUPPORT TRANSCRIPTS"

echo "This script will run STEP 19 of the DENTAL pipeline:"
echo "=== ${step} ==="; echo ""
echo "Description: Check read support per transcript and filter out unsupported"

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
#-->19 FILTER OUT LOW SUPPORT TRANSCRIPTS
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

#   19 FILTER OUT LOW SUPPORT TRANSCRIPTS
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
trinstats="${trinscripts}/util/TrinityStats.pl"
trinalign="${trinscripts}/util/align_and_estimate_abundance.pl"
trinfilter="${trinscripts}/util/filter_low_expr_transcripts.pl"
trinfilterold="/software/trinityrnaseq/trinityrnaseq-2.2.0/util/filter_fasta_by_rsem_values.pl"
trincdhit="$trinscripts/util/misc/filter_similar_seqs_expr_and_strand_aware.pl"
supertranscripts="${trinscripts}/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py"

# Modules:
bowtie2=bowtie2/2.2.9
cdhit=cd-hit/4.6.8-2017-1208
rsem=RSEM/1.3.1
samtools=samtools/1.9
seqtk=seqtk/20191028

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
mkdir ${MainDirectory}/Filter

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
memreq=64

# Memory required (in MB):
memreqmb=`echo "${memreq}*1000" | bc`

# CPUs required (integer):
cpureq=16

# Memory per thread:
mempt=`echo "${memreq}/${cpureq}"|bc`

# FPKM filter cutoffs:
HfpkmcutS=1; HfpkmcutU=5
LfpkmcutS=0.5; LfpkmcutU=3

# CDhit cutoffs:
identity=98; coverage=90

############################################
##       CREATE PBS FILES FOR FILTER      ##
############################################
# PBS files for this step are named Filter_id_PROJECT.pbs
# Cluster job names for this step are named Filter_id_PROJECT

# List all fastq files and create correspondent PBS files:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do cp ${MainDirectory}/Header.pbs Filter_${id}_${project}.pbs; done; done

# Populate PBS files (name of job):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -N Filter_${id}_${project}" >> Filter_${id}_${project}.pbs; done; done

# Populate PBS files (number of threads):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -r n" >> Filter_${id}_${project}.pbs; done; done

# Populate PBS files (time and resources for running the job):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -l mem=${memreq}GB,walltime=320:00:00,ncpus=${cpureq}" >> Filter_${id}_${project}.pbs; done; done

# Populate PBS files (send message when job is finished or aborted):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -m ae" >> Filter_${id}_${project}.pbs; done; done

# Populate PBS files (e-mail address for correspondence):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "#PBS -M $email" >> Filter_${id}_${project}.pbs; done; done

# Populate PBS files (empty line):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Filter_${id}_${project}.pbs; done; done


# Populate PBS files (empty line):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $seqtk" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $bowtie2" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $samtools" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $rsem" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "module load $cdhit" >> Filter_${id}_${project}.pbs; done; done

# Populate PBS files (empty line):
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Filter_${id}_${project}.pbs; done; done


############################################
##         COMMAND LINE FOR FILTER        ##
############################################
# Populate PBS files (actual command line):

# Filtered Assemblies:

# Create directories:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "mkdir ${MainDirectory}/FPKMcutFiltered_${id}_${project}" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "mkdir ${MainDirectory}/FPKMcutFiltered_${id}_${project}/Spliced ${MainDirectory}/FPKMcutFiltered_${id}_${project}/Unspliced" >> Filter_${id}_${project}.pbs; done; done

# Enter directory:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd ${MainDirectory}/FPKMcutFiltered_${id}_${project}" >> Filter_${id}_${project}.pbs; done; done

# Make inputs from alignment:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "grep -v "^@" ${MainDirectory}/GmapFiltered_${id}_${project}/Transrate_${id}_${project}.spliced.sam | cut -f1 >> ${MainDirectory}/GmapFiltered_${id}_${project}/Transrate_${id}_${project}.spliced.ids" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "grep -v "^@" ${MainDirectory}/GmapFiltered_${id}_${project}/Transrate_${id}_${project}.unspliced.sam | cut -f1 >> ${MainDirectory}/GmapFiltered_${id}_${project}/Transrate_${id}_${project}.unspliced.ids" >> Filter_${id}_${project}.pbs; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "seqtk subseq ${MainDirectory}/Transrate_${id}_${project}/Trinity_${id}_${project}.Trinity/good.Trinity_${id}_${project}.Trinity.fasta ${MainDirectory}/GmapFiltered_${id}_${project}/Transrate_${id}_${project}.spliced.ids >> ${MainDirectory}/FPKMcutFiltered_${id}_${project}/Transrate_${id}_${project}.spliced.fasta" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "seqtk subseq ${MainDirectory}/Transrate_${id}_${project}/Trinity_${id}_${project}.Trinity/good.Trinity_${id}_${project}.Trinity.fasta ${MainDirectory}/GmapFiltered_${id}_${project}/Transrate_${id}_${project}.unspliced.ids >> ${MainDirectory}/FPKMcutFiltered_${id}_${project}/Transrate_${id}_${project}.unspliced.fasta" >> Filter_${id}_${project}.pbs; done; done

# Filter with CDhit based on sequence similarity:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd-hit-est -o ${MainDirectory}/FPKMcutFiltered_${id}_${project}/CDhit${identity}_${id}_${project}.spliced.fasta -c 0.${identity} -i ${MainDirectory}/FPKMcutFiltered_${id}_${project}/Transrate_${id}_${project}.spliced.fasta -g 1 -r 0 -aL 0.9 -aS 0.${coverage} -p 1 -d 0 -b 3 -T ${cpureq} -M ${memreqmb}" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd-hit-est -o ${MainDirectory}/FPKMcutFiltered_${id}_${project}/CDhit${identity}_${id}_${project}.unspliced.fasta -c 0.${identity} -i ${MainDirectory}/FPKMcutFiltered_${id}_${project}/Transrate_${id}_${project}.unspliced.fasta -g 1 -r 0 -aL 0.9 -aS 0.${coverage} -p 1 -d 0 -b 3 -T ${cpureq} -M ${memreqmb}" >> Filter_${id}_${project}.pbs; done; done

# Automatically align with Bowtie2, count reads with RSEM and filter by FPKM:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${trinalign} --transcripts ${MainDirectory}/FPKMcutFiltered_${id}_${project}/CDhit${identity}_${id}_${project}.spliced.fasta --seqType fq --left ${MainDirectory}/Trinity_${id}_${project}/insilico_read_normalization/left.norm.fq.paired.fq --right ${MainDirectory}/Trinity_${id}_${project}/insilico_read_normalization/right.norm.fq.paired.fq --est_method RSEM --coordsort_bam --output_dir ${MainDirectory}/FPKMcutFiltered_${id}_${project}/Spliced --aln_method bowtie2 --SS_lib_type RF --thread_count ${cpureq} --prep_reference --trinity_mode --bowtie2_RSEM" '"'"-q --no-mixed --no-discordant --gbar 1000 --np 0 --end-to-end -k 200 -p ${cpureq}"'"' >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${trinalign} --transcripts ${MainDirectory}/FPKMcutFiltered_${id}_${project}/CDhit${identity}_${id}_${project}.unspliced.fasta --seqType fq --left ${MainDirectory}/Trinity_${id}_${project}/insilico_read_normalization/left.norm.fq.paired.fq --right ${MainDirectory}/Trinity_${id}_${project}/insilico_read_normalization/right.norm.fq.paired.fq --est_method RSEM --coordsort_bam --output_dir ${MainDirectory}/FPKMcutFiltered_${id}_${project}/Unspliced --aln_method bowtie2 --SS_lib_type RF --thread_count ${cpureq} --prep_reference --trinity_mode --bowtie2_RSEM" '"'"-q --no-mixed --no-discordant --gbar 1000 --np 0 --end-to-end -k 200 -p ${cpureq}"'"' >> Filter_${id}_${project}.pbs; done; done

# Filter by FPKM based on RSEM results:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd ${MainDirectory}/FPKMcutFiltered_${id}_${project}/Spliced; ${trinfilterold} --rsem_output=${MainDirectory}/FPKMcutFiltered_${id}_${project}/Spliced/RSEM.isoforms.results --fasta ${MainDirectory}/FPKMcutFiltered_${id}_${project}/CDhit${identity}_${id}_${project}.spliced.fasta --output ${MainDirectory}/FPKMcutFiltered_${id}_${project}/FPKM_${id}_${project}.spliced.fasta --fpkm_cutoff=${HfpkmcutS}" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd ${MainDirectory}/FPKMcutFiltered_${id}_${project}/Unspliced; ${trinfilterold} --rsem_output=${MainDirectory}/FPKMcutFiltered_${id}_${project}/Unspliced/RSEM.isoforms.results --fasta ${MainDirectory}/FPKMcutFiltered_${id}_${project}/CDhit${identity}_${id}_${project}.unspliced.fasta --output ${MainDirectory}/FPKMcutFiltered_${id}_${project}/FPKM_${id}_${project}.unspliced.fasta --fpkm_cutoff=${HfpkmcutU}" >> Filter_${id}_${project}.pbs; done; done

# Filter by Low FPKM based on RSEM results:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd ${MainDirectory}/FPKMcutFiltered_${id}_${project}/Spliced; ${trinfilterold} --rsem_output=${MainDirectory}/FPKMcutFiltered_${id}_${project}/Spliced/RSEM.isoforms.results --fasta ${MainDirectory}/FPKMcutFiltered_${id}_${project}/CDhit${identity}_${id}_${project}.spliced.fasta --output ${MainDirectory}/FPKMcutFiltered_${id}_${project}/LowFPKM_${id}_${project}.spliced.fasta --fpkm_cutoff=${LfpkmcutS}" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd ${MainDirectory}/FPKMcutFiltered_${id}_${project}/Unspliced; ${trinfilterold} --rsem_output=${MainDirectory}/FPKMcutFiltered_${id}_${project}/Unspliced/RSEM.isoforms.results --fasta ${MainDirectory}/FPKMcutFiltered_${id}_${project}/CDhit${identity}_${id}_${project}.unspliced.fasta --output ${MainDirectory}/FPKMcutFiltered_${id}_${project}/LowFPKM_${id}_${project}.unspliced.fasta --fpkm_cutoff=${LfpkmcutU}" >> Filter_${id}_${project}.pbs; done; done


# Separate:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "" >> Filter_${id}_${project}.pbs; done; done



# Raw Assemblies:

# Create and enter directory:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "mkdir ${MainDirectory}/FPKMcutRaw_${id}_${project}" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "mkdir ${MainDirectory}/FPKMcutRaw_${id}_${project}/Spliced ${MainDirectory}/FPKMcutRaw_${id}_${project}/Unspliced" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd ${MainDirectory}/FPKMcutRaw_${id}_${project}" >> Filter_${id}_${project}.pbs; done; done

# Make inputs from alignment:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "grep -v "^@" ${MainDirectory}/GmapRaw_${id}_${project}/Trinity_${id}_${project}.spliced.sam | cut -f1 >> ${MainDirectory}/GmapRaw_${id}_${project}/Trinity_${id}_${project}.spliced.ids" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "grep -v "^@" ${MainDirectory}/GmapRaw_${id}_${project}/Trinity_${id}_${project}.unspliced.sam | cut -f1 >> ${MainDirectory}/GmapRaw_${id}_${project}/Trinity_${id}_${project}.unspliced.ids" >> Filter_${id}_${project}.pbs; done; done

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "seqtk subseq ${MainDirectory}/Trinity/Trinity_${id}_${project}.Trinity.fasta ${MainDirectory}/GmapRaw_${id}_${project}/Trinity_${id}_${project}.spliced.ids >> ${MainDirectory}/FPKMcutRaw_${id}_${project}/Trinity_${id}_${project}.spliced.fasta" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "seqtk subseq ${MainDirectory}/Trinity/Trinity_${id}_${project}.Trinity.fasta ${MainDirectory}/GmapRaw_${id}_${project}/Trinity_${id}_${project}.unspliced.ids >> ${MainDirectory}/FPKMcutRaw_${id}_${project}/Trinity_${id}_${project}.unspliced.fasta" >> Filter_${id}_${project}.pbs; done; done

# Filter with CDhit based on sequence similarity:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd-hit-est -o ${MainDirectory}/FPKMcutRaw_${id}_${project}/CDhit${identity}_${id}_${project}.spliced.fasta -c 0.${identity} -i ${MainDirectory}/FPKMcutRaw_${id}_${project}/Trinity_${id}_${project}.spliced.fasta -g 1 -r 0 -aL 0.9 -aS 0.${coverage} -p 1 -d 0 -b 3 -T ${cpureq} -M ${memreqmb}" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd-hit-est -o ${MainDirectory}/FPKMcutRaw_${id}_${project}/CDhit${identity}_${id}_${project}.unspliced.fasta -c 0.${identity} -i ${MainDirectory}/FPKMcutRaw_${id}_${project}/Trinity_${id}_${project}.unspliced.fasta -g 1 -r 0 -aL 0.9 -aS 0.${coverage} -p 1 -d 0 -b 3 -T ${cpureq} -M ${memreqmb}" >> Filter_${id}_${project}.pbs; done; done


# Automatically align with Bowtie2, count reads with RSEM and filter by FPKM:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${trinalign} --transcripts ${MainDirectory}/FPKMcutRaw_${id}_${project}/CDhit${identity}_${id}_${project}.spliced.fasta --seqType fq --left ${MainDirectory}/Trinity_${id}_${project}/insilico_read_normalization/left.norm.fq.paired.fq --right ${MainDirectory}/Trinity_${id}_${project}/insilico_read_normalization/right.norm.fq.paired.fq --est_method RSEM --coordsort_bam --output_dir ${MainDirectory}/FPKMcutRaw_${id}_${project}/Spliced --aln_method bowtie2 --SS_lib_type RF --thread_count ${cpureq} --prep_reference --trinity_mode --bowtie2_RSEM" '"'"-q --no-mixed --no-discordant --gbar 1000 --np 0 --end-to-end -k 200 --quiet -p ${cpureq}"'"' >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${trinalign} --transcripts ${MainDirectory}/FPKMcutRaw_${id}_${project}/CDhit${identity}_${id}_${project}.unspliced.fasta --seqType fq --left ${MainDirectory}/Trinity_${id}_${project}/insilico_read_normalization/left.norm.fq.paired.fq --right ${MainDirectory}/Trinity_${id}_${project}/insilico_read_normalization/right.norm.fq.paired.fq --est_method RSEM --coordsort_bam --output_dir ${MainDirectory}/FPKMcutRaw_${id}_${project}/Unspliced --aln_method bowtie2 --SS_lib_type RF --thread_count ${cpureq} --prep_reference --trinity_mode --bowtie2_RSEM" '"'"-q --no-mixed --no-discordant --gbar 1000 --np 0 --end-to-end -k 200 --quiet -p ${cpureq}"'"' >> Filter_${id}_${project}.pbs; done; done

# Filter by FPKM based on RSEM results:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd ${MainDirectory}/FPKMcutRaw_${id}_${project}/Spliced; ${trinfilterold} --rsem_output=${MainDirectory}/FPKMcutRaw_${id}_${project}/Spliced/RSEM.isoforms.results --fasta ${MainDirectory}/FPKMcutRaw_${id}_${project}/CDhit${identity}_${id}_${project}.spliced.fasta --output ${MainDirectory}/FPKMcutRaw_${id}_${project}/FPKM_${id}_${project}.spliced.fasta --fpkm_cutoff=${HfpkmcutS}" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd ${MainDirectory}/FPKMcutRaw_${id}_${project}/Unspliced; ${trinfilterold} --rsem_output=${MainDirectory}/FPKMcutRaw_${id}_${project}/Unspliced/RSEM.isoforms.results --fasta ${MainDirectory}/FPKMcutRaw_${id}_${project}/CDhit${identity}_${id}_${project}.unspliced.fasta --output ${MainDirectory}/FPKMcutRaw_${id}_${project}/FPKM_${id}_${project}.unspliced.fasta --fpkm_cutoff=${HfpkmcutU}" >> Filter_${id}_${project}.pbs; done; done

# Filter by Low FPKM based on RSEM results:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd ${MainDirectory}/FPKMcutRaw_${id}_${project}/Spliced; ${trinfilterold} --rsem_output=${MainDirectory}/FPKMcutRaw_${id}_${project}/Spliced/RSEM.isoforms.results --fasta ${MainDirectory}/FPKMcutRaw_${id}_${project}/CDhit${identity}_${id}_${project}.spliced.fasta --output ${MainDirectory}/FPKMcutRaw_${id}_${project}/LowFPKM_${id}_${project}.spliced.fasta --fpkm_cutoff=${LfpkmcutS}" >> Filter_${id}_${project}.pbs; done; done
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "cd ${MainDirectory}/FPKMcutRaw_${id}_${project}/Unspliced; ${trinfilterold} --rsem_output=${MainDirectory}/FPKMcutRaw_${id}_${project}/Unspliced/RSEM.isoforms.results --fasta ${MainDirectory}/FPKMcutRaw_${id}_${project}/CDhit${identity}_${id}_${project}.unspliced.fasta --output ${MainDirectory}/FPKMcutRaw_${id}_${project}/LowFPKM_${id}_${project}.unspliced.fasta --fpkm_cutoff=${LfpkmcutU}" >> Filter_${id}_${project}.pbs; done; done


############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (header):
echo "" >> $MainDirectory/Pipeline_Log
echo "===> Read support" >> $MainDirectory/Pipeline_Log

# Populate Log file (empty line):
echo "" >> $MainDirectory/Pipeline_Log

# Populate Log file (date):
echo "Starting at:" >> $MainDirectory/Pipeline_Log
date >> $MainDirectory/Pipeline_Log
echo "" >> $MainDirectory/Pipeline_Log
echo "-- command lines --" >> $MainDirectory/Pipeline_Log

# Populate Log file (command lines used):
tail Filter_*.pbs >> $MainDirectory/Pipeline_Log

############################################
##           SUBMIT FILTER JOBS           ##
############################################
# Populate Log file (header):
echo "Job Submission IDs:" >> ${MainDirectory}/PBSsub/SubJobs_Filter

# Submit jobs and populate Log file (submission ID numbers):
ls Filter_*.pbs | while read file; do echo "$file" >> ${MainDirectory}/PBSsub/SubJobs_Filter ; qsub "$file"; done >> ${MainDirectory}/PBSsub/SubJobs_Filter

# Inform status of jobs after submission:
echo "Your ${step} jobs are currently running:"
qstat | grep "$user" | grep Filter
qstat | grep -c Filter

############################################
##                   END                  ##
############################################

