# Create output directory:
mkdir AllOutMGFR

# Create 500 input files for MGFR:
seq 1 500 | while read n; do echo "Basal Her2 LumA LumB Normal" | tr ' ' '\n' | while read subtype; do grep $subtype FullSampleIDs_withType.sorted | sort -R | head -n3 >> ThreeReplicates_$n; done; done

seq 1 500 | while read n; do cut -f1 ThreeReplicates_$n | while read id; do grep -n $id FullSampleIDs >> ColNumbers_$n; done; done

seq 1 500 | while read n; do sort -nb -k1,1 ColNumbers_$n | cut -d':' -f2- | while read id; do grep $id ThreeReplicates_$n >> InputMGFR_$n; done; done

seq 1 500 | while read n; do cols=`cut -d':' -f1 ColNumbers_$n | tr '\n' ','`; cut -d',' -f${cols}1 RSEM_AllTCGACount.csv >> RSEM_AllTCGACount_$n.csv; done

seq 1 500 | while read n; do cols=`cut -d':' -f1 ColNumbers_$n | tr '\n' ','`; cut -d',' -f${cols}1 RSEM_AllTCGACount_NoAnnotated.csv >> RSEM_AllTCGACount_NoAnnotated_$n.csv; done

seq 1 500 | while read n; do cols=`cut -d':' -f1 ColNumbers_$n | tr '\n' ','`; cut -d',' -f${cols}1 RSEM_AllTCGACount_AnnotatedOnly.csv >> RSEM_AllTCGACount_AnnotatedOnly_$n.csv; done


# Load R:
module load R/4.0.2 

# Run MGFR for all inputs:
seq 1 500 | while read n; do Rscript MGFR.R RSEM_AllTCGACount_${n}.csv InputMGFR_${n} OutMGFR_${n}; mv OutMGFR_${n}* AllOutMGFR; echo "===== $n ====="; done

# Run MGFR for all annotated gene inputs:
seq 1 500 | while read n; do Rscript MGFR.R RSEM_AllTCGACount_AnnotatedOnly_${n}.csv InputMGFR_${n} OutMGFR_AnnotatedOnly_${n}; mv OutMGFR_AnnotatedOnly_${n}* AllOutMGFR; echo "===== $n ====="; done

# Run MGFR for all NB-lncRNA inputs:
seq 1 500 | while read n; do Rscript MGFR.R RSEM_AllTCGACount_NoAnnotated_${n}.csv InputMGFR_${n} OutMGFR_NoAnnotated_${n}; mv OutMGFR_NoAnnotated_${n}* AllOutMGFR; echo "===== $n ====="; done

mkdir MGFR500Rounds; mv ColNumbers_* MGFR500Rounds; mv InputMGFR_* MGFR500Rounds; mv ThreeReplicates_* MGFR500Rounds; mv AllOutMGFR AllOutMGFR500Rounds


