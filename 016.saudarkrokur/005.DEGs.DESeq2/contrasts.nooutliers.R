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
library(ramify)

#
# 0. user-defined variables
#
setwd("~/scratch/")
kallisto_dir = "/Users/adrian/research/016.saudarkrokur/results/kallisto/kallisto.dme.100"
results_dir = '/Users/adrian/research/016.saudarkrokur/results/deseq2'

#
# 1. generate gene to transcript mapping
#
listEnsembl()
listEnsembl(version=113)
#ensembl = useEnsembl(biomart="ensembl", verbose=TRUE)
ensembl = useEnsembl(biomart="ensembl", verbose=TRUE, mirror='asia')
head(listDatasets(ensembl)) # dmelanogaster_gene_ensembl
mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                        dataset="dmelanogaster_gene_ensembl",
                        #host = 'https://www.ensembl.org',
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
conditions = c(rep('h1M8', 3), rep('h2F14', 3), rep('h2M8', 3), rep('ko', 3), rep('wt', 3))

metadata = data.frame(labels)
metadata$condition = conditions
metadata$path = paths
View(metadata)

#
# 3. contrasts
#
read_threshold = 20
effect_size_threshold = log2(2)
tpm_threshold = 2

contrasts = list()
contrasts[[1]] = c('wt', 'ko')
contrasts[[2]] = c('h1M8', 'ko')
contrasts[[3]] = c('h2F14', 'wt')
contrasts[[4]] = c('h2M8', 'wt')

contrast_maker <- function(contrast){
  
  label = paste(contrast[1], contrast[2], sep='_')
  message(label)
  
  rule = (metadata$condition == contrast[1]) | (metadata$condition == contrast[2])
  working_metadata = metadata[rule, ]
  
  outlier_samples = c('h2F14_3', 'h2M8_2')
  for (forbidden in outlier_samples){
    if (forbidden %in% working_metadata$labels){
      print('found')
      index = which(working_metadata$labels == forbidden)
      working_metadata = working_metadata[-index, ]
    }
  }
  
  print(dim(working_metadata))
  print(working_metadata)
  
  txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)
  dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~condition) 
  dds$condition = relevel(dds$condition, contrast[2])
  
  # keep features with at least 20 counts median difference.
  cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
  if (dim(working_metadata)[1] == 6){}
  
  
  
  a = counts(dds)[ , 1:3]
  b = counts(dds)[ , 4:6]
  c = rowMedians(a) - rowMedians(b)
  keep = abs(c) >= read_threshold
  dds = dds[keep, ]
  cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
  
  # keep features with at least a max median expression of 1 TPM.
  cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
  subset = txi$abundance[names(dds), ]
  a = rowMedians(subset[ , 1:3])
  b = rowMedians(subset[ , 4:6])
  c = pmax(a, b)
  keep = c >= tpm_threshold
  dds = dds[keep, ]
  cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
  
  # run test
  dds = DESeq(dds, test="LRT", reduced=~1)
  
  # handle tests
  res = results(dds, parallel=TRUE, alpha=0.05) 
  filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
  sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
  anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
  cat(blue(paste('contrast DEGs:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
  
  # add annotation and expression values
  subset = txi$abundance[names(dds), ]
  a = rowMedians(subset[ , 1:3])
  b = rowMedians(subset[ , 4:6])
  sorted_filtred_results[paste('expression', contrast[1], sep='')] = a
  sorted_filtred_results[paste('expression', contrast[2], sep='')] = b
  df_new = t2g[t2g$ensembl_gene_id %in% rownames(sorted_filtred_results), ]
  sorted_filtred_results['description'] = df_new$description
  sorted_filtred_results['external_gene_name'] = df_new$external_gene_name
  sorted_filtred_results['gene_biotype'] = df_new$gene_biotype
  sorted_filtred_results['entrezgene_id'] = df_new$entrezgene_id
  
  write.table(sorted_filtred_results, file=paste(results_dir, '/effect_', label, '.tsv', sep=''), quote=FALSE, sep='\t')
  write.table(anti_results, file=paste(results_dir, '/effect_', label, '.anti.tsv', sep=''), quote=FALSE, sep='\t')
  
  #               
  # volcano
  #
  plotting_x = sorted_filtred_results$log2FoldChange
  y = sorted_filtred_results$padj
  epsilon = min(y[y !=0])
  plotting_y = -log10(y + epsilon) 
  z = log10(rowMedians(txi$abundance[rownames(sorted_filtred_results), ]) + 1)
  df = data.frame(plotting_x=clip(plotting_x, .min=-6, .max=6), plotting_y=clip(plotting_y, .min=0, .max=20), plotting_z=clip(z, .min=0, .max=3))
  reds = df[df$plotting_x > 0, ]
  blues = df[df$plotting_x < 0, ]
  
  plotting_x = anti_results$log2FoldChange
  plotting_y = -log10(anti_results$padj)
  blacks = data.frame(plotting_x=clip(plotting_x, .min=-6, .max=6), plotting_y=clip(plotting_y, .min=0, .max=20))
  
  ggplot() + 
    geom_point(data=reds, aes(x=plotting_x, y=plotting_y, color=plotting_z), , size=3, shape=19, alpha=2/3, stroke=0) + 
    geom_point(data=blues, aes(plotting_x, plotting_y, color=plotting_z), size=3, shape=19, alpha=2/3, stroke=0) +
    geom_point(data=blacks, aes(plotting_x, plotting_y), color = "black", size=1, shape=19, alpha=0.2, stroke=0) +
    labs(x=expression('Expression [log'[2]~'FC]'), y=expression('Significance [log'[10]~'adjusted P]'), color=expression('Expression [log'[10]~'TPM]'), title=label) + 
    theme_linedraw() +
    geom_segment(aes(x=-1, xend=-1, y=-log10(0.05), yend=20), linetype=2) +
    geom_segment(aes(x=1, xend=1, y=-log10(0.05), yend=20), linetype=2) +
    geom_segment(aes(x=-6, xend=-1, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
    geom_segment(aes(x=1, xend=6, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
    xlim(-6.2, 6) +
    scale_color_viridis_c(option = "cividis") 
  ggsave(paste(label, '.png', sep=''))

  message('...')
}

for (contrast in contrasts){
  contrast_maker(contrast)
}