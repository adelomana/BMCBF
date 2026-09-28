rm(list = ls())

#
# -1. install libraries
# 
#if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("rhdf5")
# 
# BiocManager::install("DESeq2")
# BiocManager::install("tximport")

#
# 0. load libraries
#
library(DESeq2)
library(tximport)
library(stringr)
library(rhdf5)

#
# 1. user-defined variables
#
kallisto_dir = "/Users/adrian/research/095/quant"
results_dir = kallisto_dir

#
# 1. generate gene to transcript mapping
#
df = read.csv('/Users/adrian/research/095/ref/t2g.txt', sep='\t', header=FALSE)
t2g = df
dim(t2g)

#
# 2. define metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
labels = sapply(strsplit(dirnames, split='/',fixed=TRUE), function(x) (x[7]))
labels = str_remove(labels, '_processed')
print(labels)
paths = file.path(dirnames, 'abundance.h5')

metadata = data.frame(labels)
metadata$path = paths
View(metadata)

#
# 3. read files
#
txi = tximport(metadata$path, type="kallisto", tx2gene=t2g)

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
store = paste(results_dir, '/DESeq2_TPM_values.tsv', sep='')
write.table(tpm, file=store, quote=FALSE, sep='\t', col.names=NA)
