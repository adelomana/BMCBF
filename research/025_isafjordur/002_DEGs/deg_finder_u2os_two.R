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
library(stringr)
library(ramify) # this is for clip, for the volcano

#
# 0. user-defined variables
#
setwd("~/scratch/")
kallisto_dir = "/Users/adrian/research/bmcbf/025_isafjordur/results/000_quantification"
results_dir = '/Users/adrian/research/bmcbf/025_isafjordur/results/002_DEGs'

#
# 1. generate gene to transcript mapping
#
df = read.csv('/Users/adrian/software/kallisto/human_index_standard/annotation.tsv', sep='\t')
t2g = df[, 2:3]
dim(t2g)

#
# 2. define metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
dirnames = dirnames[grep('U2', dirnames)]
dirnames = dirnames[grep('UT', dirnames)]
paths = file.path(dirnames, 'kallisto_output_c/abundance.h5')
paths = paths[c(1, 4, 5, 6)]

labels = sapply(strsplit(paths, split='/',fixed=TRUE), function(x) (x[9]))
labels = str_remove(labels, '_processed')

genotypes = rep(c('wt', 'ko'), 2)

metadata = data.frame(labels)
metadata$path = paths
metadata$genotype = genotypes
View(metadata)

#
# 3. contrasts
#
threshold = 20
effect_size_threshold = log2(2)
tpm_threshold = 1

#
# 3.1. contrast 
#
txi = tximport(metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)
dds = DESeqDataSetFromTximport(txi, colData=metadata, design=~genotype) 
dds$genotype = relevel(dds$genotype, ref="wt")

# keep features with at least 20 counts median difference
cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
a = counts(dds)[ , 1:2]
b = counts(dds)[ , 3:4]
c = rowMedians(a) - rowMedians(b)
keep = abs(c) >= threshold
dds = dds[keep, ]
cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)

# keep features with at least a max median expression of 1 TPM.
cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
subset = txi$abundance[names(dds), ]
a = rowMedians(subset[ , 1:2])
b = rowMedians(subset[ , 3:4])
c = pmax(a, b)
keep = c >= tpm_threshold
dds = dds[keep, ]
cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) # it does not seem to affect  https://www.biostars.org/p/209118/ 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast wt vs control:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, file=paste(results_dir, '/effect_ko_vs_wt.u2os.two.for.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, file=paste(results_dir, '/effect_ko_vs_wt.u2os.two.anti.tsv', sep=''), quote=FALSE, sep='\t')
write.table(res, file=paste(results_dir, '/effect_ko_vs_wt.u2os.two.full.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('genotype')) + ggtitle('effect ko vs wt | two')

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
  labs(x=expression('Expression [log'[2]~'FC]'), y=expression('Significance [log'[10]~'adjusted P]'), color=expression('Expression [log'[10]~'TPM]'), title='KO vs WT | U2OS two') + 
  theme_linedraw() +
  geom_segment(aes(x=-1, xend=-1, y=-log10(0.05), yend=20), linetype=2) +
  geom_segment(aes(x=1, xend=1, y=-log10(0.05), yend=20), linetype=2) +
  geom_segment(aes(x=-6, xend=-1, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  geom_segment(aes(x=1, xend=6, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  xlim(-6, 6) +
  scale_color_viridis_c(option = "cividis") +
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=24))
ggsave('u2os.png')
dev.off()
