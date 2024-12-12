rm(list = ls())

#
# -1. install libraries
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.20")
# 
# BiocManager::install("biomaRt")
# BiocManager::install("tximport")
# BiocManager::install("DESeq2")
# BiocManager::install("rhdf5")

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
kallisto_dir = "/Users/adrian/research/016.saudarkrokur/results/kallisto/kallisto.dme.100"
results_dir = '/Users/adrian/research/016.saudarkrokur/results/deseq2'

#
# 1. generate gene to transcript mapping
#
listEnsembl()
listEnsembl(version=113)
ensembl = useEnsembl(biomart="ensembl")
head(listDatasets(ensembl)) # dmelanogaster_gene_ensembl
mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                        dataset="dmelanogaster_gene_ensembl",
                        host = 'https://www.ensembl.org',
                        verbose = TRUE)
# attributes = listAttributes(mart)
# hgnc_symbol gives less than external_gene_name
working_attributes = c('ensembl_transcript_id', 
                      'ensembl_gene_id', 
                      'external_gene_name', 
                      'entrezgene_id',
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
condition = c(rep('h1M8', 3), rep('h2F14', 3), rep('h2M8', 3), rep('KO', 3), rep('WT', 3))

metadata = data.frame(labels)
metadata$condition = condition
metadata$path = paths
View(metadata)

#
# 3. read files
#
txi = tximport(metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

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

store = paste(results_dir, '/annotation.tsv', sep='')
write.table(t2g, file=store, quote=FALSE, sep='\t', col.names=NA)
