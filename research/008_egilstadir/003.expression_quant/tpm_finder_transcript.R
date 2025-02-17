rm(list = ls())

#
# -1. install libraries
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# BiocManager::install("tximport")

#
# 0. load libraries
#
library(DESeq2)
library(tximport)
library(biomaRt)

#
# 1. user-defined variables
#
setwd("~/scratch/")
kallisto_dir = "/Users/adrian/research/egilsstadir/results/kallisto/kallisto.100"
results_dir = '/Users/adrian/research/egilsstadir/results/deseq2'

#
# 2. define metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
paths = file.path(dirnames, 'abundance.h5')
labels = sapply(strsplit(paths, split='/',fixed=TRUE), function(x) (x[9]))
metadata = data.frame(labels)
metadata$path = paths
View(metadata)

#
# 3. read files
#
txi = tximport(metadata$path, type="kallisto", txOut=TRUE)

#
# 4. find abundance
#
tpm = txi$abundance
colnames(tpm) = metadata$labels
dim(tpm)
View(tpm)

#
# 5. store
#
store = paste(results_dir, '/DESeq2_TPM_values_transcript.tsv', sep='')
write.table(tpm, file=store, quote=FALSE, sep='\t', col.names=NA)
