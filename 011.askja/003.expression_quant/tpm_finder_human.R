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
kallisto_dir = "/Users/adrian/research/011.askja/results/kallisto/kallisto.100"
results_dir = '/Users/adrian/research/011.askja/results/deseq2'

#
# 1. generate gene to transcript mapping
#
mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                        dataset="hsapiens_gene_ensembl",
                        host = 'https://www.ensembl.org',
                        verbose = TRUE)
# attributes = listAttributes(mart)
# hgnc_symbol gives less than external_gene_name
working_attributes = c('ensembl_transcript_id', 
                      'ensembl_gene_id', 
                      'external_gene_name',
                      'gene_biotype',
                      'description')
t2g = biomaRt::getBM(attributes=working_attributes, 
                     mart=mart,
                     verbose=TRUE)
dim(t2g)
View(t2g)

#
# 2. define metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
paths = file.path(dirnames, 'abundance.h5')
labels = sapply(strsplit(paths, split='/',fixed=TRUE), function(x) (x[9]))
replicates = rep(c('A', 'B', 'C'), 8)
timepoints = rep(c(rep('D120', 3), rep('D240', 3), rep('D060', 3), rep('control', 3)), 2)
species = c(rep('dog', 12), rep('human', 12))

metadata = data.frame(labels)
metadata$replicate = replicates
metadata$timepoint = timepoints
metadata$species = species
metadata$path = paths
View(metadata)

working_metadata = metadata[metadata$species == 'human', ]
working_metadata = working_metadata[order(working_metadata$timepoint), ]

#
# 3. read files
#
txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

#
# 4. find abundance
#
tpm = txi$abundance
colnames(tpm) = working_metadata$labels
dim(tpm)
View(tpm)

#
# 5. store
#
store = paste(results_dir, '/DESeq2_TPM_values.human.tsv', sep='')
write.table(tpm, file=store, quote=FALSE, sep='\t', col.names=NA)

store = paste(results_dir, '/annotation.human.tsv', sep='')
write.table(t2g, file=store, quote=FALSE, sep='\t', col.names=NA)
