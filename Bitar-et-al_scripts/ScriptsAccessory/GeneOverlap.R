#!/usr/bin/env Rscript

args<-commandArgs(TRUE)
library("GeneOverlap")

A <- read.table(args[1])
B <- read.table(args[2])

go.obj <- newGeneOverlap(A$V1, B$V1, genome.size = args[3])
go.obj <- testGeneOverlap(go.obj)
sink(args[4])
print(go.obj)
sink()

