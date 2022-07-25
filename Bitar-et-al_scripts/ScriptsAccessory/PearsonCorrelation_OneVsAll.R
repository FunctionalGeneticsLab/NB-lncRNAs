					#!/usr/bin/env Rscript
# module load R/4.0.2

library(tidyverse)

# args[1] = Input file
# args[2] = Column name
# args[3] = Output file (root)

args<-commandArgs(TRUE)

data <- read.delim(args[1], header=T, sep="\t")

pvvector <- c()
corvector <- c()
namevector <- c()
colname <- as.character(args[2])

for (name in colnames(data)[!colnames(data) %in% "FileID"]) {
	pvcor <- cor.test(data[[colname]], data[[name]], na.action = "na.exclude", method="pearson")$p.value
	pvvector <- c(pvvector, pvcor)
	cor <- cor(data[[colname]], use = "complete.obs", data[[name]], method="pearson")
	corvector <- c(corvector, cor)
	namevector <- c(namevector, name)
}

df <- data.frame(Names=as.character(namevector), Cov=as.numeric(corvector), Pval=as.numeric(pvvector)) 

write.table(df, file=sprintf("PearsonCorrelations_%s_%s.table",args[3],args[2]),sep="\t", quote=FALSE)

