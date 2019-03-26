#! /usr/bin/env Rscript
library(ape)
args = commandArgs(trailingOnly=TRUE)

seqs = read.FASTA(args[1])
cat(c(
    dist.dna(seqs, model='raw')[1], 
    dist.dna(seqs, model='JC69')[1], 
    dist.dna(seqs, model='TN93')[1],
    dist.dna(seqs, model='K80')[1],
    dist.dna(seqs, model='T92')[1]
), sep='\t', fill=TRUE)

