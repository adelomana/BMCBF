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
kallisto_dir = "/Users/adrian/research/egilsstadir/results/kallisto/kallisto.100"
results_dir = '/Users/adrian/research/egilsstadir/results/deseq2'

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
View(t2g)

#
# 2. define metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
paths = file.path(dirnames, 'abundance.h5')
labels = sapply(strsplit(paths, split='/',fixed=TRUE), function(x) (x[9]))
condition = c(rep('202high', 3), rep('202low', 3), rep('KO', 3), rep('WT', 3))

metadata = data.frame(labels)
metadata$condition = condition
metadata$path = paths
View(metadata)

#
# 3. contrasts
#
threshold = 10
effect_size_threshold = log2(2)

#
# 3.1. contrast WT vs KO
#
rule = (metadata$condition == 'KO') | (metadata$condition == 'WT')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~condition) 
dds$time = relevel(dds$condition, ref="KO")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast WT vs KO:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/effect_WT_vs_KO.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/effect_WT_vs_KO.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('condition')) + ggtitle('effect WT vs KO')
#ggsave(file.path(results_dir, 'effect_WT_vs_KO.png'))

#
# 3.2. contrasts high vs KO
#
rule = (metadata$condition == '202high') | (metadata$condition == 'KO')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~condition) 
dds$time = relevel(dds$condition, ref="KO")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast 202high vs KO:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/effect_202high_vs_KO.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/effect_202high_vs_KO.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('condition')) + ggtitle('effect 202high vs KO')
#ggsave(file.path(results_dir, 'effect_202high_vs_WT.png'))

#
# 3.3. contrasts low vs KO
#
rule = (metadata$condition == '202low') | (metadata$condition == 'KO')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~condition) 
dds$time = relevel(dds$condition, ref="KO")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast 202low vs KO:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/effect_202low_vs_KO.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/effect_202low_vs_KO.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('condition')) + ggtitle('effect 202low vs KO')
#ggsave(file.path(results_dir, 'effect_202low_vs_KO.png'))

#
# 3.3. contrasts high vs low
#
rule = (metadata$condition == '202high') | (metadata$condition == '202low')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~condition) 
dds$time = relevel(dds$condition, ref="202low")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast 202high vs 202low:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/effect_202high_vs_202low.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/effect_202high_vs_202low.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('condition')) + ggtitle('effect high vs low')
#ggsave(file.path(results_dir, 'effect_202high_vs_202low.png'))

