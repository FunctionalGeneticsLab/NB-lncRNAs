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
cd $MainDirectory

# TransRate score (0-1) captures how confident you can be in what was assembled, as well as how complete the assembly is.
echo "" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
echo "########## OVERALL SCORE ##########" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${id}_${project}" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; head -n1 ${MainDirectory}/Transrate_${id}_${project}/assemblies.csv | tr ',' '\n' >> tempfileA; tail -n1 ${MainDirectory}/Transrate_${id}_${project}/assemblies.csv | tr ',' '\n' >> tempfileB; paste tempfileA tempfileB | grep "^score" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; rm -rf tempfileA tempfileB; echo "####################" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; done; done

# The expression-weighted assembly score is more generous to assemblies with many low-expressed contigs that are badly assembled.
# Only works in version 1.0.3 (not version 1.0.1)
#echo "" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
#echo "########## WEIGHTED SCORE ##########" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
#echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${id}_${project}" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; head -n1 ${MainDirectory}/Transrate_${id}_${project}/assemblies.csv | tr ',' '\n' >> tempfileA; tail -n1 ${MainDirectory}/Transrate_${id}_${project}/assemblies.csv | tr ',' '\n' >> tempfileB; paste tempfileA tempfileB | grep "weighted" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; rm -rf tempfileA tempfileB; echo "####################" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; done; done

# TransRate can automatically assess the best cutoff and use contig scores to filter out bad contigs from an assembly, keeping only well-assembled ones to maximise (optimise) the assembly score.
echo "" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
echo "########## OPTIMAL SCORE ##########" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${id}_${project}" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; head -n1 ${MainDirectory}/Transrate_${id}_${project}/assemblies.csv | tr ',' '\n' >> tempfileA; tail -n1 ${MainDirectory}/Transrate_${id}_${project}/assemblies.csv | tr ',' '\n' >> tempfileB; paste tempfileA tempfileB | grep "optimal_score" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; rm -rf tempfileA tempfileB; echo "####################" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; done; done

# The contig score can be thought of as measure of whether the contig is an accurate, complete, non-redundant representation of a transcript that was present in the sequenced sample.
# Contig scores can be used to filter out bad contigs from an assembly, leaving you with only the well-assembled ones.
echo "" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
echo "########## CONTIG SCORE ##########" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${id}_${project}" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; cat "${MainDirectory}/Transrate_${id}_${project}/Trinity_${id}_${project}.Trinity/contigs.csv" | tr ',' '\t' | cut -f9 | "$StatScript" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; done; done

#'Good' mappings are those where both members of the pair are aligned, in the correct orientation, on the same contig and without overlapping either end of the contig.
echo "" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
echo "########## REFERENCE COVERAGE ##########" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${id}_${project}" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; head -n1 ${MainDirectory}/Transrate_${id}_${project}/assemblies.csv | tr ',' '\n' >> tempfileA; tail -n1 ${MainDirectory}/Transrate_${id}_${project}/assemblies.csv | tr ',' '\n' >> tempfileB; paste tempfileA tempfileB | grep "reference_coverage" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; rm -rf tempfileA tempfileB; echo "####################" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; done; done

# The percentage of reads that map to the reference:
echo "" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
echo "########## PERCENTAGE MAPPED ##########" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${id}_${project}" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; head -n1 ${MainDirectory}/Transrate_${id}_${project}/assemblies.csv | tr ',' '\n' >> tempfileA; tail -n1 ${MainDirectory}/Transrate_${id}_${project}/assemblies.csv | tr ',' '\n' >> tempfileB; paste tempfileA tempfileB | grep "p_fragments_mapped" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; rm -rf tempfileA tempfileB; echo "####################" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; done; done

#'Good' mappings are those where both members of the pair are aligned, in the correct orientation, on the same contig and without overlapping either end of the contig.
echo "" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
echo "########## GOOD MAPPINGS % ##########" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${id}_${project}" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; head -n1 ${MainDirectory}/Transrate_${id}_${project}/assemblies.csv | tr ',' '\n' >> tempfileA; tail -n1 ${MainDirectory}/Transrate_${id}_${project}/assemblies.csv | tr ',' '\n' >> tempfileB; paste tempfileA tempfileB | grep "p_good_mapping" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; rm -rf tempfileA tempfileB; echo "####################" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; done; done

#'Good' mappings are those where both members of the pair are aligned, in the correct orientation, on the same contig and without overlapping either end of the contig.
echo "" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
echo "########## MEAN ORF % ##########" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${id}_${project}" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; head -n1 ${MainDirectory}/Transrate_${id}_${project}/assemblies.csv | tr ',' '\n' >> tempfileA; tail -n1 ${MainDirectory}/Transrate_${id}_${project}/assemblies.csv | tr ',' '\n' >> tempfileB; paste tempfileA tempfileB | grep "mean_orf_percent" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; rm -rf tempfileA tempfileB; echo "####################" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; done; done

echo "" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
echo "########## % OF BASES UNCOVERED ##########" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate;
echo "A,B" | tr ',' '\n' | while read id; do cut -d'_' -f1 ListOfFiles | sort | uniq | while read project; do echo "${id}_${project}" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; head -n1 ${MainDirectory}/Transrate_${id}_${project}/assemblies.csv | tr ',' '\n' >> tempfileA; tail -n1 ${MainDirectory}/Transrate_${id}_${project}/assemblies.csv | tr ',' '\n' >> tempfileB; paste tempfileA tempfileB | grep "p_bases_uncovered" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; rm -rf tempfileA tempfileB; echo "####################" >> ${MainDirectory}/StepSummaries/SummaryOfTransRate; done; done

# Or simply report everything:
#cut -d'_' -f1 ListOfFiles | uniq | while read folder; do echo "" >> ${MainDirectory}/Trinity/fullSummaryOfTransRate; echo "$folder" >> ${MainDirectory}/Trinity/fullSummaryOfTransRate ; cat ${MainDirectory}/TransRate_${folder}.o* >> ${MainDirectory}/Trinity/fullSummaryOfTransRate; done

############################################
##            WRITE TO LOG FILE           ##
############################################
# Populate Log file (TransRate results):
echo "** TransRate Results **" >> ${MainDirectory}/Pipeline_Log

cut -d'_' -f1 ListOfFiles | uniq | while read folder; do echo "" >> ${MainDirectory}/Pipeline_Log; echo "$folder" >> ${MainDirectory}/Pipeline_Log ; cat ${MainDirectory}/Busco_${folder}/short_summary_Busco_${folder}.txt >> ${MainDirectory}/Pipeline_Log; done

############################################
##          MOVE FILES TO OUTPUT          ##
############################################
mv TransRate_*.pbs ${MainDirectory}/PBSin
mv *TransRate_*.e* ${MainDirectory}/PBSout
mv *TransRate_*.o* ${MainDirectory}/PBSout
cd $MainDirectory
############################################
##             STEP  FINISHED             ##
############################################


