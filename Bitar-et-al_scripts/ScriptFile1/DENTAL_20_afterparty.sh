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
##         CHECK STEP CONCLUSION!         ##
############################################

############################################
##         SET MAIN WORK DIRECTORY        ##
############################################
# Define directory:

if [ "$#" -gt 0 ]; then MainDirectory=$1; fi

if [ "$#" -eq 0 ]; then MainDirectory=`pwd`; fi

echo ""; echo "WARNING: The main directory for this run was set to ${MainDirectory}."; echo ""

############################################
##              MAKE SUMMARY              ##
############################################

cd ${MainDirectory}

# Create auxiliary awk function (once only):
#echo "{" >> ${MainDirectory}/count.awk ; echo '  print length($0);' >> ${MainDirectory}/count.awk; echo "}" >> ${MainDirectory}/count.awk

# In-house Statistics script:
StatScript="/working/lab_julietF/mainaB/Statistics.sh"

# Retrieve length and exon number for each SuperTranscript:

### RAW files:
#echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do grep "^>" ${MainDirectory}/SuperTranscripts/rawSuper_${id}_${project}.fasta | tr -d '>' >> ${MainDirectory}/SuperTranscripts/TEMP1_${id}_${project}; grep -v "^>" ${MainDirectory}/SuperTranscripts/rawSuper_${id}_${project}.fasta | awk -f count.awk >> ${MainDirectory}/SuperTranscripts/TEMP2_${id}_${project}; cut -f1 ${MainDirectory}/SuperTranscripts/rawSuper_${id}_${project}.gtf | grep "TRINITY" | uniq -c | cut -c1-8 | tr -d ' ' >> ${MainDirectory}/SuperTranscripts/TEMP3_${id}_${project}; echo -e "SuperTranscriptID\tLength\tExonNumber" >> ${MainDirectory}/SuperTranscripts/rawSuper_${id}_${project}.metrics; paste ${MainDirectory}/SuperTranscripts/TEMP1_${id}_${project} ${MainDirectory}/SuperTranscripts/TEMP2_${id}_${project} ${MainDirectory}/SuperTranscripts/TEMP3_${id}_${project} >> ${MainDirectory}/SuperTranscripts/rawSuper_${id}_${project}.metrics; rm -rf ${MainDirectory}/SuperTranscripts/TEMP1_${id}_${project} ${MainDirectory}/SuperTranscripts/TEMP2_${id}_${project} ${MainDirectory}/SuperTranscripts/TEMP3_${id}_${project}; done; done

#echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do grep "^>" ${MainDirectory}/SuperTranscripts/rawLowSuper_${id}_${project}.fasta | tr -d '>' >> ${MainDirectory}/SuperTranscripts/TEMP1_${id}_${project}; grep -v "^>" ${MainDirectory}/SuperTranscripts/rawLowSuper_${id}_${project}.fasta | awk -f count.awk >> ${MainDirectory}/SuperTranscripts/TEMP2_${id}_${project}; cut -f1 ${MainDirectory}/SuperTranscripts/rawLowSuper_${id}_${project}.gtf | grep "TRINITY" | uniq -c | cut -c1-8 | tr -d ' ' >> ${MainDirectory}/SuperTranscripts/TEMP3_${id}_${project}; echo -e "SuperTranscriptID\tLength\tExonNumber" >> ${MainDirectory}/SuperTranscripts/rawLowSuper_${id}_${project}.metrics; paste ${MainDirectory}/SuperTranscripts/TEMP1_${id}_${project} ${MainDirectory}/SuperTranscripts/TEMP2_${id}_${project} ${MainDirectory}/SuperTranscripts/TEMP3_${id}_${project} >> ${MainDirectory}/SuperTranscripts/rawLowSuper_${id}_${project}.metrics; rm -rf ${MainDirectory}/SuperTranscripts/TEMP1_${id}_${project} ${MainDirectory}/SuperTranscripts/TEMP2_${id}_${project} ${MainDirectory}/SuperTranscripts/TEMP3_${id}_${project}; done; done


### FILTERED files:
#echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do grep "^>" ${MainDirectory}/SuperTranscripts/filteredSuper_${id}_${project}.fasta | tr -d '>' >> ${MainDirectory}/SuperTranscripts/TEMP1_${id}_${project}; grep -v "^>" ${MainDirectory}/SuperTranscripts/filteredSuper_${id}_${project}.fasta | awk -f count.awk >> ${MainDirectory}/SuperTranscripts/TEMP2_${id}_${project}; cut -f1 ${MainDirectory}/SuperTranscripts/filteredSuper_${id}_${project}.gtf | grep "TRINITY" | uniq -c | cut -c1-8 | tr -d ' ' >> ${MainDirectory}/SuperTranscripts/TEMP3_${id}_${project}; echo -e "SuperTranscriptID\tLength\tExonNumber" >> ${MainDirectory}/SuperTranscripts/filteredSuper_${id}_${project}.metrics; paste ${MainDirectory}/SuperTranscripts/TEMP1_${id}_${project} ${MainDirectory}/SuperTranscripts/TEMP2_${id}_${project} ${MainDirectory}/SuperTranscripts/TEMP3_${id}_${project} >> ${MainDirectory}/SuperTranscripts/filteredSuper_${id}_${project}.metrics; rm -rf ${MainDirectory}/SuperTranscripts/TEMP1_${id}_${project} ${MainDirectory}/SuperTranscripts/TEMP2_${id}_${project} ${MainDirectory}/SuperTranscripts/TEMP3_${id}_${project}; done; done

#echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do grep "^>" ${MainDirectory}/SuperTranscripts/filteredLowSuper_${id}_${project}.fasta | tr -d '>' >> ${MainDirectory}/SuperTranscripts/TEMP1_${id}_${project}; grep -v "^>" ${MainDirectory}/SuperTranscripts/filteredLowSuper_${id}_${project}.fasta | awk -f count.awk >> ${MainDirectory}/SuperTranscripts/TEMP2_${id}_${project}; cut -f1 ${MainDirectory}/SuperTranscripts/filteredLowSuper_${id}_${project}.gtf | grep "TRINITY" | uniq -c | cut -c1-8 | tr -d ' ' >> ${MainDirectory}/SuperTranscripts/TEMP3_${id}_${project}; echo -e "SuperTranscriptID\tLength\tExonNumber" >> ${MainDirectory}/SuperTranscripts/filteredLowSuper_${id}_${project}.metrics; paste ${MainDirectory}/SuperTranscripts/TEMP1_${id}_${project} ${MainDirectory}/SuperTranscripts/TEMP2_${id}_${project} ${MainDirectory}/SuperTranscripts/TEMP3_${id}_${project} >> ${MainDirectory}/SuperTranscripts/filteredLowSuper_${id}_${project}.metrics; rm -rf ${MainDirectory}/SuperTranscripts/TEMP1_${id}_${project} ${MainDirectory}/SuperTranscripts/TEMP2_${id}_${project} ${MainDirectory}/SuperTranscripts/TEMP3_${id}_${project}; done; done


# Final summary:
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "### ${id} ${project} ###" >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; echo "# Raw Assembly #" >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; echo "Length Summary:" | tr '\n' ' ' >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; grep -v "SuperTranscriptID" ${MainDirectory}/SuperTranscripts/rawSuper_${id}_${project}.metrics | cut -f2 | "$StatScript" >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; echo "ExonNumber Summary:" | tr '\n' ' ' >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; grep -v "SuperTranscriptID" ${MainDirectory}/SuperTranscripts/rawSuper_${id}_${project}.metrics | cut -f3 | "$StatScript" >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; echo "# Raw-Low Assembly #" >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; echo "Length Summary:" | tr '\n' ' ' >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; grep -v "SuperTranscriptID" ${MainDirectory}/SuperTranscripts/rawLowSuper_${id}_${project}.metrics | cut -f2 | "$StatScript" >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; echo "ExonNumber Summary:" | tr '\n' ' ' >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; grep -v "SuperTranscriptID" ${MainDirectory}/SuperTranscripts/rawLowSuper_${id}_${project}.metrics | cut -f3 | "$StatScript" >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; echo "# Filtered Assembly #" >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; echo "Length Summary:" | tr '\n' ' ' >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; grep -v "SuperTranscriptID" ${MainDirectory}/SuperTranscripts/filteredSuper_${id}_${project}.metrics | cut -f2 | "$StatScript" >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; echo "ExonNumber Summary:" | tr '\n' ' ' >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; grep -v "SuperTranscriptID" ${MainDirectory}/SuperTranscripts/filteredSuper_${id}_${project}.metrics | cut -f3 | "$StatScript" >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; echo "# Filtered-Low Assembly #" >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; echo "Length Summary:" | tr '\n' ' ' >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; grep -v "SuperTranscriptID" ${MainDirectory}/SuperTranscripts/filteredLowSuper_${id}_${project}.metrics | cut -f2 | "$StatScript" >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; echo "ExonNumber Summary:" | tr '\n' ' ' >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; grep -v "SuperTranscriptID" ${MainDirectory}/SuperTranscripts/filteredLowSuper_${id}_${project}.metrics | cut -f3 | "$StatScript" >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; echo "" >> ${MainDirectory}/StepSummaries/SummaryOfSuperTranscripts; done; done

############################################
##          MOVE FILES TO OUTPUT          ##
############################################
mv Supertranscripts_*.pbs ${MainDirectory}/PBSin
mv *Supertranscripts_*.e* ${MainDirectory}/PBSout
mv *Supertranscripts_*.o* ${MainDirectory}/PBSout
cd $MainDirectory


############################################
##              MAKE SUMMARY              ##
############################################
# Summarize counts through filtering steps:

echo -e "Assembly\tTotalCounts\tLowFPKMcut\tHighFPKMcut\tTransRate\tfilteredLowFPKMcut\tfilteredHighFPKMcut" >> ${MainDirectory}/StepSummaries/SummaryOfTranscriptCounts; echo -e "Assembly\tTotalCounts\tLowFPKMcut\tHighFPKMcut\tTransRate\tfilteredLowFPKMcut\tfilteredHighFPKMcut" >> ${MainDirectory}/StepSummaries/SummaryOfGeneCounts

echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${id}_${project}"; alltrans=`grep -c "^>" ${MainDirectory}/Trinity/Trinity_${id}_${project}.Trinity.fasta`; allgenes=`grep "^>" ${MainDirectory}/Trinity/Trinity_${id}_${project}.Trinity.fasta | cut -d'_' -f1-4 | sort | uniq | wc -l`; tratetrans=`grep -c "^>" ${MainDirectory}/Transrate_${id}_${project}/Trinity_${id}_${project}.Trinity/good.Trinity_${id}_${project}.Trinity.fasta`; trategenes=`grep "^>" ${MainDirectory}/Transrate_${id}_${project}/Trinity_${id}_${project}.Trinity/good.Trinity_${id}_${project}.Trinity.fasta | cut -d'_' -f1-4 | sort | uniq | wc -l`; lffpkmtrans=`grep -c "^>" ${MainDirectory}/SuperTranscripts/filteredLowFPKM_${id}_${project}.fasta`; lffpkmgenes=`grep "^>" ${MainDirectory}/SuperTranscripts/filteredLowFPKM_${id}_${project}.fasta | cut -d'_' -f1-4 | sort | uniq | wc -l`; hffpkmtrans=`grep -c "^>" ${MainDirectory}/SuperTranscripts/filteredFPKM_${id}_${project}.fasta`; hffpkmgenes=`grep "^>" ${MainDirectory}/SuperTranscripts/filteredFPKM_${id}_${project}.fasta | cut -d'_' -f1-4 | sort | uniq | wc -l`; lrfpkmtrans=`grep -c "^>" ${MainDirectory}/SuperTranscripts/rawLowFPKM_${id}_${project}.fasta`; lrfpkmgenes=`grep "^>" ${MainDirectory}/SuperTranscripts/rawLowFPKM_${id}_${project}.fasta | cut -d'_' -f1-4 | sort | uniq | wc -l`; hrfpkmtrans=`grep -c "^>" ${MainDirectory}/SuperTranscripts/rawFPKM_${id}_${project}.fasta`; hrfpkmgenes=`grep "^>" ${MainDirectory}/SuperTranscripts/rawFPKM_${id}_${project}.fasta | cut -d'_' -f1-4 | sort | uniq | wc -l`; echo -e "${project}_${id}\t${alltrans}\t${lrfpkmtrans}\t${hrfpkmtrans}\t${tratetrans}\t${lffpkmtrans}\t${hffpkmtrans}" >> ${MainDirectory}/StepSummaries/SummaryOfTranscriptCounts; echo -e "${project}_${id}\t${allgenes}\t${lrfpkmgenes}\t${hrfpkmgenes}\t${trategenes}\t${lffpkmgenes}\t${hffpkmgenes}" >> ${MainDirectory}/StepSummaries/SummaryOfGeneCounts; done; done


############################################
##             STEP  FINISHED             ##
############################################




