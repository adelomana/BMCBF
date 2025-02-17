rm(list = ls())

#
# -1. install libraries
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")

#
# 0. load libraries
#
library(DESeq2)
library(tximport)
library(biomaRt)
library(BiocParallel)
library(crayon) 
library(ggplot2)

#
# 0. user-defined variables
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
working_attributes = c('ensembl_transcript_id', 
                       'ensembl_gene_id', 
                       'external_gene_name',
                       'gene_biotype',
                       'description')
t2g = biomaRt::getBM(attributes=working_attributes, 
                     mart=mart,
                     verbose=TRUE)
dim(t2g)

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

working_metadata = metadata[metadata$species == 'human', ]
working_metadata = working_metadata[order(working_metadata$timepoint), ]
metadata = working_metadata
View(metadata)

#
# 3. contrasts
#
threshold = 20
effect_size_threshold = log2(2)

#
# 3.1. contrast D60 vs control
#
rule = (metadata$timepoint == 'D060') | (metadata$timepoint == 'control')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~timepoint) 
dds$time = relevel(dds$timepoint, ref="control")

# keep features with at least 20 counts median difference
cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
a = counts(dds)[ , 1:3]
b = counts(dds)[ , 4:6]
c = rowMedians(a) - rowMedians(b)
keep = abs(c) >= threshold
dds = dds[keep, ]
cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast D060 vs control:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, file=paste(results_dir, '/effect_D060_vs_control.human.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, file=paste(results_dir, '/effect_D060_vs_control.anti.human.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('timepoint')) + ggtitle('effect D060 vs control')

#
# 3.2. contrast D120 vs control
#
rule = (metadata$timepoint == 'D120') | (metadata$timepoint == 'control')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~timepoint) 
dds$time = relevel(dds$timepoint, ref="control")

# keep features with at least 20 counts median difference
cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
a = counts(dds)[ , 1:3]
b = counts(dds)[ , 4:6]
c = rowMedians(a) - rowMedians(b)
keep = abs(c) >= threshold
dds = dds[keep, ]
cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast D120 vs control:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, file=paste(results_dir, '/effect_D120_vs_control.human.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, file=paste(results_dir, '/effect_D120_vs_control.anti.human.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('timepoint')) + ggtitle('effect D120 vs control')

#
# 3.3. contrast D240 vs control
#
rule = (metadata$timepoint == 'D240') | (metadata$timepoint == 'control')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~timepoint) 
dds$time = relevel(dds$timepoint, ref="control")

# keep features with at least 20 counts median difference
cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
a = counts(dds)[ , 1:3]
b = counts(dds)[ , 4:6]
c = rowMedians(a) - rowMedians(b)
keep = abs(c) >= threshold
dds = dds[keep, ]
cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast D240 vs control:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, file=paste(results_dir, '/effect_D240_vs_control.human.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, file=paste(results_dir, '/effect_D240_vs_control.anti.human.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('timepoint')) + ggtitle('effect D240 vs control')

