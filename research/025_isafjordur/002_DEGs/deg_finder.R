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
kallisto_dir = '/Users/adrian/research/019_fossvogur/001_profiles/'
results_dir = '/Users/adrian/research/011.askja/results/deseq2'

#
# 1. generate gene to transcript mapping
#
mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                        dataset="hsapiens_gene_ensembl",
                        host = 'https://oct2022.archive.ensembl.org', # because of 108
                        verbose = TRUE)
# attributes = listAttributes(mart)
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
tpm_threshold = 1

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

# keep features with at least a max median expression of 1 TPM.
cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
subset = txi$abundance[names(dds), ]
a = rowMedians(subset[ , 1:3])
b = rowMedians(subset[ , 4:6])
c = pmax(a, b)
keep = c >= tpm_threshold
dds = dds[keep, ]
cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.1) # alpha 0.1 or 0.05?  
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast D060 vs control:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, file=paste(results_dir, '/effect_D060_vs_control.human.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, file=paste(results_dir, '/effect_D060_vs_control.anti.human.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('timepoint')) + ggtitle('effect D060 vs control')














# run test
dds = DESeq(dds, test="LRT", reduced=~1)

# handle tests
res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast DEGs:', dim(filtred_results)[1], sep=' ')), fill=TRUE)

# add annotation and expression values
subset = txi$abundance[rownames(sorted_filtred_results), ]
a = rowMedians(subset[ , 1:3])
b = rowMedians(subset[ , 4:6])

sorted_filtred_results[paste('expression', contrast[1], sep='')] = a
sorted_filtred_results[paste('expression', contrast[2], sep='')] = b
tempo = t2g[t2g[ , 'ensembl_gene_id'] %in% rownames(sorted_filtred_results), ]
df_new = tempo[!duplicated(tempo$ensembl_gene_id), ]
no = match(rownames(sorted_filtred_results), df_new$ensembl_gene_id)
sorted_filtred_results['ensembl'] = df_new$ensembl_gene_id[no]
sorted_filtred_results['description'] = df_new$description[no]
sorted_filtred_results['gene_biotype'] = df_new$gene_biotype[no]
sorted_filtred_results['entrezgene_id'] = df_new$entrezgene_id[no]

write.table(sorted_filtred_results, file=paste(results_dir, '/effect_', label, '.for.tsv', sep=''), quote=FALSE, sep='\t')
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