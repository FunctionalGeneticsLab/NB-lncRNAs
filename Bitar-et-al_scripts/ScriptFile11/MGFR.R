

# Use arguments:
args<-commandArgs(TRUE)

library(MGFR)
data <- read.csv(args[1], sep=",", header=TRUE, row.names=1)

CV <- read.csv(args[2], sep="", header=FALSE)
cv <- CV$V2

M <- getMarkerGenes.rnaseq(data, class.vec=cv, samples2compare="all", annotate=FALSE, score.cutoff=0.75)

outfilebasal <- paste(args[3],"Basal",sep="_")
outfileluma <- paste(args[3],"LumA",sep="_")
outfilelumb <- paste(args[3],"LumB",sep="_")
outfilenormal <- paste(args[3],"Normal",sep="_")
outfileher2 <- paste(args[3],"Her2",sep="_")

write.table(M$Basal,file = outfilebasal, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(M$LumA,file = outfileluma, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(M$LumB,file = outfilelumb, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(M$Normal,file = outfilenormal, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(M$Her2,file = outfileher2, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)


