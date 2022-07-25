#!/usr/bin/env Rscript
# module load R/4.0.2

library(genefu)

args<-commandArgs(TRUE)

annotation <- read.delim(as.matrix(args[1], header=T, sep='\t'))

rnaseq <- read.delim(as.matrix(args[2], header=T, sep='\t'))

# Prepare TPM input:
n <- rnaseq$probe
rnaseq_t <- as.data.frame(t(rnaseq[,-1]))
colnames(rnaseq_t) <- n

# PAM50:
data(pam50)
pam50_predictions_all <- molecular.subtyping(sbt.model = "pam50", data = rnaseq_t, annot = annotation, do.mapping = T)
write.table(pam50_predictions_all$subtype, file=sprintf("Pam50Classification_%s",args[3]), sep='\t', quote=F)
write.table(pam50_predictions_all$subtype.proba, file=sprintf("Pam50ClassificationProbabilities_%s",args[3]), sep='\t', quote=F)
write.table(pam50_predictions_all$subtype.crisp, file=sprintf("Pam50ClassificationCrisp_%s",args[3]), sep='\t', quote=F)

# Claudin Low:
data(claudinLowData)
claudin_predictions_all <- molecular.subtyping(sbt.model = "claudinLow", data = rnaseq_t, annot = annotation, do.mapping = T)
write.table(claudin_predictions_all$subtype, file=sprintf("ClaudinClassification_%s",args[3]), sep='\t', quote=F)
write.table(claudin_predictions_all$subtype.proba, file=sprintf("ClaudinClassificationProbabilities_%s",args[3]), sep='\t', quote=F)
write.table(claudin_predictions_all$subtype.crisp, file=sprintf("ClaudinClassificationCrisp_%s",args[3]), quote=F)




