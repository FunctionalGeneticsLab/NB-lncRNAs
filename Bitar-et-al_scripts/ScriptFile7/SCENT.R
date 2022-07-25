#Source:
#https://rdrr.io/github/aet21/SCENT/f/vignettes/SCENT.Rmd

# Input File (for example):
#cd /working/lab_julietF/mainaB/Project_NormalBreast/scRNAseq/Individuals1to3/BowtieRsem_GencodeAll/Stemness

module load R/4.0.2

library(SCENT);
data(net13Jun12);
data(net17Jan16);
print(dim(net13Jun12.m));
# [1] 8434 8434
print(dim(net17Jan16.m));
# [1] 11751 11751

pheno.v <- unlist(read.table("CellClusters/Sorted_CellClusters_Input", sep="\t", header=FALSE))

netlncp.c  <- as.matrix(read.table("/working/lab_julietF/mainaB/StemnessNetworks/LncRNAs2Pathways/SCENTinput/LncPathNetwork.csv", sep=",", header=TRUE, row.names=1))
netlncp.m <- apply(netlncp.c, 2, as.numeric)
rownames(netlncp.m) <- rownames(netlncp.c)
print(dim(netlncp.m)); typeof(netlncp.m)


entrez.m <- read.csv("RSEM_AllClustersCount_inbuiltSCENTinput.csv", header=T, row.names=1, sep=",")
lncp.m <- read.csv("RSEM_AllClustersCount_lncpathSCENTinput.csv", header=T, row.names=1,sep=",")
lncp.m[is.na(lncp.m)] <- 0.00

logentrez.m <- log2(entrez.m+1.1)
loglncp.m <- log2(lncp.m+1.1)

intlogentrez.net13 <- DoIntegPPI(exp.m = logentrez.m, ppiA.m = net13Jun12.m)
intrawentrez.net13 <- DoIntegPPI(exp.m = entrez.m, ppiA.m = net13Jun12.m)

intlogentrez.net17 <- DoIntegPPI(exp.m = logentrez.m, ppiA.m = net17Jan16.m)
intrawentrez.net17 <- DoIntegPPI(exp.m = entrez.m, ppiA.m = net17Jan16.m)

intloglncp.lncp <- DoIntegPPI(exp.m = loglncp.m, ppiA.m = netlncp.m)
intrawlncp.lncp <- DoIntegPPI(exp.m = lncp.m, ppiA.m = netlncp.m)

#intloglncp.lncp <- as.matrix(intloglncp.lncp.df)

str(intlogentrez.net13)
str(intrawentrez.net13)

str(intlogentrez.net17)
str(intrawentrez.net17)

str(intloglncp.lncp)
str(intrawlncp.lncp)


intloglncp.lncp.sr <- CompSRana(intloglncp.lncp, local = FALSE, mc.cores = 30)
intrawlncp.lncp.sr <- CompSRana(intrawlncp.lncp, local = FALSE, mc.cores = 30)

intlogentrez.net13.sr <- CompSRana(intlogentrez.net13, local = FALSE, mc.cores = 30)
intrawentrez.net13.sr <- CompSRana(intrawentrez.net13, local = FALSE, mc.cores = 30)

intlogentrez.net17.sr <- CompSRana(intrawentrez.net17, local = FALSE, mc.cores = 30)
intrawentrez.net17.sr <- CompSRana(intrawentrez.net17, local = FALSE, mc.cores = 30)


intlogentrez.net13.v <- intlogentrez.net13.sr$SR
intrawentrez.net13.v <- intrawentrez.net13.sr$SR

intlogentrez.net17.v <- intlogentrez.net17.sr$SR
intrawentrez.net17.v <- intrawentrez.net17.sr$SR

intloglncp.lncp.v <- intloglncp.lncp.sr$SR
intrawlncp.lncp.v <- intrawlncp.lncp.sr$SR


write.table(intlogentrez.net13.v, file="SR_Net13_log2_SR", sep="\n", quote=F)
write.table(intrawentrez.net13.v, file="SR_Net13_raw_SR", sep="\n", quote=F)

write.table(intlogentrez.net17.v, file="SR_Net17_log2_SR", sep="\n", quote=F)
write.table(intrawentrez.net17.v, file="SR_Net17_raw_SR", sep="\n", quote=F)

write.table(intloglncp.lncp.v, file="SR_LncP_log2_SR", sep="\n", quote=F)
write.table(intrawlncp.lncp.v, file="SR_LncP_raw_SR", sep="\n", quote=F)


jpeg(file="SR_LncP_log2.jpeg")
boxplot(intloglncp.lncp.v ~ pheno.v, main = "SR potency estimates", xlab = "Cell Type", ylab = "SR")
dev.off()

jpeg(file="SR_LncP_raw.jpeg")
boxplot(intrawlncp.lncp.v ~ pheno.v, main = "SR potency estimates", xlab = "Cell Type", ylab = "SR")
dev.off()

jpeg(file="SR_Net13_log2.jpeg")
boxplot(intlogentrez.net13.v ~ pheno.v, main = "SR potency estimates", xlab = "Cell Type", ylab = "SR")
dev.off()

jpeg(file="SR_Net13_raw.jpeg")
boxplot(intrawentrez.net13.v ~ pheno.v, main = "SR potency estimates", xlab = "Cell Type", ylab = "SR")
dev.off()

jpeg(file="SR_Net17_log2.jpeg")
boxplot(intlogentrez.net17.v ~ pheno.v, main = "SR potency estimates", xlab = "Cell Type", ylab = "SR")
dev.off()

jpeg(file="SR_Net17_raw.jpeg")
boxplot(intrawentrez.net17.v ~ pheno.v, main = "SR potency estimates", xlab = "Cell Type", ylab = "SR")
dev.off()


intlogentrez.net13.pot <- InferPotencyStates(potest.v=intlogentrez.net13.v, pheno.v = pheno.v)
intrawentrez.net13.pot <- InferPotencyStates(potest.v=intrawentrez.net13.v, pheno.v = pheno.v)

intlogentrez.net17.pot <- InferPotencyStates(potest.v=intlogentrez.net17.v, pheno.v = pheno.v)
intrawentrez.net17.pot <- InferPotencyStates(potest.v=intrawentrez.net17.v, pheno.v = pheno.v)

intloglncp.lncp.pot <- InferPotencyStates(potest.v=intloglncp.lncp.v, pheno.v = pheno.v)
intrawlncp.lncp.pot <- InferPotencyStates(potest.v=intrawlncp.lncp.v, pheno.v = pheno.v)



write.table(intlogentrez.net13.pot$distr, file="SR_Net13_log2_PotSum", sep="\t", quote=F)
write.table(intrawentrez.net13.pot$distr, file="SR_Net13_raw_PotSum", sep="\t", quote=F)

write.table(intlogentrez.net17.pot$distr, file="SR_Net17_log2_PotSum", sep="\t", quote=F)
write.table(intrawentrez.net17.pot$distr, file="SR_Net17_raw_PotSum", sep="\t", quote=F)

write.table(intloglncp.lncp.pot$distr, file="SR_LncP_log2_PotSum", sep="\t", quote=F)
write.table(intrawlncp.lncp.pot$distr, file="SR_LncP_raw_PotSum", sep="\t", quote=F)


intlogentrez.net13.pot$distr
intlogentrez.net17.pot$distr
intloglncp.lncp.pot$distr

