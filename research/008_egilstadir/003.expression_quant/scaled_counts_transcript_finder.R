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
raw_counts = txi$counts
colnames(raw_counts) = metadata$labels
dim(raw_counts)
View(raw_counts)

factor = colSums(txi$counts)/mean(colSums(txi$counts))
print(factor)
scaled_counts = raw_counts/factor
View(scaled_counts)

print(raw_counts[1, ])
print(factor)
print(scaled_counts[1, ])

#
# 5. store
#
store = paste(results_dir, '/DESeq2_TPM_values_transcript.counts.tsv', sep='')
write.table(scaled_counts, file=store, quote=FALSE, sep='\t', col.names=NA)
